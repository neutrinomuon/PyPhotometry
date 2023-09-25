! ###########################################################################
!     RESUME : Compute photometry in filters, using different               !
!              magnitude systems. Original routine dates back to 2003-2004. !
!                                                                           !
!              Original routines first released:                            !
!                                                                           !
!              Wed Sep 15 09:32:47 WEST 2004                                !
!                                                                           !
!     INPUT    : 01) O_lambda  -> Wavelength                                !
!                02) O_fluxes  -> Fluxes                                    !
!                03) T_lambda  -> Wavelength                                !
!                04) T_fluxes  -> Transmission curve                        !
!                05) file_out  -> In case of error                          !
!                06) T_l_Area  -> Area of filter (lambda)                   !
!                07) T_n_Area  -> Area of filter (frequency)                !
!                08) magABsys  -> Calibration for AB system                 !
!                09) magTGsys  -> Calibration for TG system                 !
!                10) standard  -> standard flux for VEGA                    !
!                11) lamb_eff  -> Effective wavelength                      !
!                12) Numb_lbd  -> Number of points in filter                !
!                13) name_fil  -> Name of bandpass/filter                   !
!                14) NOlambda  -> Number of O_lambda                        !
!                15) Nfilters  -> Number of filters                         !
!                16) Int_Type  -> Intergral type                            !
!                17) IsKeepOn  -> Variable:       0 => Run aborted          !
!                18) verbosity -> Optional variable to print and check      !
!                                                                           !
!     OUTPUT   : 01) MAG_spec -> Absolute magnitudes in various systems     !
!                                                                           !
!     EXTRA ROUTINES : EvalTransmission and IntegralALL                     !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!     Checked: Mon Sep 25 10:07:32 PM WEST 2023                             !
!                                                                           !
!     LOG: 25/09/2023: There was an error in the f2py declaration of        !
!     O_lambda, O_fluxes and NOlambda.                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvalFilters( O_lambda,O_fluxes,NOlambda,T_lambda,T_fluxes,       &
                        Ntlambda,Nfilters,MAG_spec,PhotFlux,Mcalibra,       &
                        T_l_Area,magABsys,magTGsys,standard,Numb_lbd,       &
                        IsKeepOn,Int_Type,LSun_con,verbosity )
    use ModDataType
    implicit none

    real     (kind=RP), parameter :: cl_speed=299792.458000_RP
    real     (kind=RP), parameter :: fac4Pid2=1.19649518e40_RP
    integer  (kind=IB), parameter :: Ncalibra=7

    ! Astronomy related lsun = 3.839e26 * watt = Lsun
    ! Within pyphot library the value is set == 3.82600000e33_RP! 3.8300500000e33_RP !=3.916e33
    ! And with neutrinos !3.828e+33 !3.82600000e33_RP

    real     (kind=RP), optional :: LSun_con!=3.839e33! This could be different from this version
    
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: NOlambda,Nfilters,Ntlambda!,Ncalibra
    integer  (kind=IB), optional :: Int_Type,verbosity
    integer  (kind=IB), dimension(Nfilters), intent(in) :: Numb_lbd
    integer  (kind=IB) :: IsShowOn,i_filter,j_filter,N_lambda,IntegraT

    real     (kind=RP), target, dimension(Ntlambda,Nfilters), intent(in) :: &
                                                                  T_lambda, &
                                                                  T_fluxes

    real     (kind=RP), pointer, dimension(:) :: P1lambda,P1fluxes
    real     (kind=RP), dimension(Ncalibra), intent(out) :: Mcalibra

    real     (kind=RP), dimension(Nfilters,Ncalibra), intent(out) :: MAG_spec

    real     (kind=RP), dimension(NOlambda), intent(in) :: O_lambda,O_fluxes
    
    real     (kind=RP), dimension(Nfilters), intent(in) :: T_l_Area,        &
                                                           magABsys,        &
                                                           magTGsys,        &
                                                           standard

    real     (kind=RP), dimension(Nfilters), intent(out) :: PhotFlux

    real     (kind=RP) :: aux_area,flux_T01,mag_Sgal,fluxVega,LSun

    character (len=CH), allocatable, dimension(:) :: Fnamecal
    character (len=CH) :: W1aux,W2aux,W3aux

    !f2py real     (kind=RP), intent(in)  :: O_lambda
    !f2py real     (kind=RP), intent(in)  :: O_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(O_fluxes) :: NOlambda=shape(O_fluxes,0)

    !f2py real     (kind=RP), intent(in)  :: T_lambda
    !f2py real     (kind=RP), intent(in)  :: T_fluxes
    !f2py real     (kind=RP), intent(in)  :: T_l_Area
    !f2py real     (kind=RP), intent(in)  :: magABsys
    !f2py real     (kind=RP), intent(in)  :: magTGsys
    !f2py real     (kind=RP), intent(in)  :: standard
    !f2py integer  (kind=IB), intent(in)  :: Numb_lbd
    !f2py real     (kind=RP), intent(out) :: PhotFlux
    !f2py real     (kind=RP), intent(out) :: Mcalibra
    !f2py real     (kind=RP), intent(out) :: MAG_spec
    !f2py integer  (kind=IB), intent(hide), depend(T_lambda) :: Ntlambda=shape(T_lambda,0), Nfilters=shape(T_lambda,1)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn

    !f2py  real     (kind=RP), optional :: LSun_con=3.839e33
    !f2py  integer  (kind=IB), optional :: Int_Type=2
    !f2py  integer  (kind=IB), optional :: verbosity=0

    interface
       subroutine EvalTransmission( L_lambda,S_fluxes,N_lambda,T_lambda,    &
                                    T_fluxes,Ntlambda,Int_Type,fluxtran,    &
                                    IsKeepOn,verbosity )
         use ModDataType
         implicit none
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(N_lambda), intent(in) :: L_lambda,   &
                                                                S_fluxes
         real     (kind=RP), dimension(Ntlambda), intent(in) :: T_lambda,   &
                                                                T_fluxes
         real     (kind=RP), intent(out) :: fluxtran
       end subroutine EvalTransmission
       subroutine EvalTFluxes( O_lambda,O_fluxes,Nspeclbd,T_lambda,T_fluxes,&
                               N_lambda,T_l_area,fluxtran,lamb_eff,IskeepOn,&
                               verbosity )
         use ModDataType
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: N_lambda,Nspeclbd
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(Nspeclbd), intent(in) :: O_lambda,   &
                                                                O_fluxes
         real     (kind=RP), dimension(N_lambda), intent(in) :: T_lambda,   &
                                                                T_fluxes
         real     (kind=RP), intent(out) :: lamb_eff,fluxtran
         real     (kind=RP), intent(in) :: T_l_area
       end subroutine EvalTFluxes
       real (kind=RP) function IntegralALL( SXvalues,SYvalues,lambda_i,      &
                                            lambda_f,N_lambda,IsKeepOn,      &
                                            Int_Type,verbosity )
         use ModDataType
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: N_lambda,Int_Type
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,    &
                                                                SYvalues
         real     (kind=RP), intent(in) :: lambda_i,lambda_f
       end function IntegralALL
    end interface

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if
    
    if ( present(Int_Type) ) then
       IntegraT = Int_Type
    else
       IntegraT = 2_IB
    end if

    if ( present(LSun_con) ) then
       LSun = LSun_con
    else
       LSun = 3.839e33_RP
    end if

    ! DEBUG due to f2py declaration
    ! integer  (kind=IB), intent(hide), depend(NOlambda) :: NOlambda=shape(O_fluxes,0)
    ! and it should be
    ! integer  (kind=IB), intent(hide), depend(O_fluxes) :: NOlambda=shape(O_fluxes,0)
    !
    !write(*,*) '=====>', O_lambda(1),size(O_lambda),O_lambda(size(O_lambda)),O_lambda(NOlambda),NOlambda
    !return
    
    allocate( Fnamecal(Ncalibra) )
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)')   '[EvalFilters]'
    end if
    
    do i_filter=1,Nfilters
       N_lambda = Numb_lbd(i_filter)
       if (associated (P1lambda)) nullify( P1lambda )
       if (associated (P1fluxes)) nullify( P1fluxes )
       
       P1lambda => T_lambda(1:N_lambda,i_filter)
       P1fluxes => T_fluxes(1:N_lambda,i_filter)
       
       aux_area = T_l_Area(i_filter)
       
! *** Evaluate transmission flux ********************************************
       !call EvalTFluxes( O_lambda,O_fluxes,P1lambda,P1fluxes,aux_area,     &
       !                  N_lambda,flux_T01,lamb_T01,NOlambda,IskeepOn,     &
       !                  IsShowOn )
       call EvalTransmission( O_lambda,O_fluxes,NOlambda,P1lambda,P1fluxes, &
                              N_lambda,IntegraT,flux_T01,IsKeepOn,IsShowOn )
! *** Evaluate transmission flux ********************************************

! *** Photometric flux ******************************************************
!     It assumes the units of input spectrum - f(lambda)                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !write (*,*) 'flux_T01 =====>',flux_T01,IsKeepOn
       PhotFlux(i_filter) = flux_T01 !* LSun / fac4Pid2
! *** Photometric flux ******************************************************

! *** Evaluate magnitudes ***************************************************
!     RESUME : Absolute magnitudes are evaluated if the object is           !
!              moved 10 pc away. The flux assumed is in erg/s,              !
!              therefore the constants are:                                 !
!                                                                           !
!              01 pc    = 3.08567758 x 10^18 cm, then                       !
!                                                                           !
!              10 pc    = 3.08567758 x 10^19 cm and this implies            !
!                                                                           !
!              4 Pi d^2 = 1.19649518 x 10^40 cm^2                           !
!                                                                           !
!              WARNING : The synthetic fluxes must be in erg/s/A, but       !
!              the spectrum is usually evaluated in L_sun/A.                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Modified because it could be here the problem in the log expression
       !flux_T01 = flux_T01 * (LSun / fac4Pid2)
       mag_Sgal = -2.50_RP * (log10(flux_T01)+log10(LSun)-log10(fac4Pid2))
       fluxVega = standard(i_filter)      ! Store value to compute ZERO point

! *** Compute ZERO points in filters ****************************************

! *** VEGA standard ************************************************************************ !
!     Ref.: [1] Bessel M.S., 2005 Annu. Rev. Astrophys., 43, 293                             !
!           [2] Cousins A.W.J., Jones D.H.P, 1976 Mem. R. astr. Soc, 81, 1.                  !
!           [3] Kitchin C.R., 2003 book, Astrophysical Techniques, ISBN: 0-7503-0946-6       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( fluxVega > 0.0_RP ) then
          Mcalibra(1) = +2.5_RP * log10(fluxVega)
       else
          Mcalibra(1) = -999.0_RP
       end if
       Fnamecal(1) = 'VEGA_std'

! *** VEGA proposed by Bohlin and Gilliland 2004 ******************************************* !
!     Ref.: Bohlin, R. C.; Gilliland, R. L. 2004, AJ, 128, 3053B                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( fluxVega > 0.0_RP ) then
          Mcalibra(2) = +2.5_RP * log10(fluxVega) + 0.026_RP
       else
          Mcalibra(2) = -999.0_RP
       end if
       Fnamecal(2) = 'VEGABG04'

! *** M_AB standard system (based on frequency) ******************************************** !
!     Ref.: [1] Oke, J.B. 1974, ApJS, 27, 21                                                 !
!                                                                                            !
!     Formula:  AB(nu) = -2.5 log10 f(nu) - 48.60                                            !
!                                                                                            !
!               lambda*f(lambda) = nu*f(nu) thus,                                            !
!                                                                                            !
!               f(nu) = f(lambda) * (lambda/nu) = f(lambda)*(lambda)^2 / c                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( fluxVega > 0.0_RP .AND. magABsys(i_filter) > -999.0_RP ) then 
          Mcalibra(3) = +2.5_RP * log10(fluxVega) + magABsys(i_filter)
       else
          Mcalibra(3) = -999.0_RP
       end if
       Fnamecal(3) = 'ABsystem'

! *** M_TG standard system - Thuan and Gunn (BD+17d4708) *********************************** !
!     Ref.: [1] Oke, J. B., and Gunn, J. E. 1983, ApJ, 266, 713                              !
!           [2] Schild, R. 1984, ApJ, 286, 450                                               !
!           [3] Schneider, D. P., Gunn, J. E., and Hoessel J. G. 1983, ApJ, 264, 337         !
!           [4] Thuan, T. X., and Gunn, J. E. 1976, PASP, 88, 543                            !
!           [5] Wade, R. A., Hoessel, J. G., Elias, J. H., Huchra, J. P. 1979, PASP, 91, 35  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( magTGsys(i_filter) > -999.0_RP ) then
          Mcalibra(4) = +2.5_RP * log10(fluxVega) + magTGsys(i_filter)
       else
          Mcalibra(4) = -999.0_RP
       endif
       Fnamecal(4) = 'TGsystem'
       
! --- WFPC2 system ------------------------------------------------------------------------- !
!     Ref.: [1] Stone, R.P.S. 1996, ApJS, 107, 423                                           !
!                                                                                            !
!     Formula: STMAG(lambda) = -2.5 log10 f(lambda) - 21.1                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Mcalibra(5) = -21.100_RP
       Fnamecal(5) = 'WFPC2sys'

! --- FOCA at 2000 ----------------------------------------------------------
       Mcalibra(6) = -21.175_RP
       Fnamecal(6) = 'FOCA_sys'

! --- Without any calibration -----------------------------------------------
       Mcalibra(7) = +0.0000_RP
       Fnamecal(7) = 'NO_calib'

! *** Compute ZERO points in filters ****************************************
       if ( IsShowOn == 1_IB ) then
          !write (*,'(4x,a,a)') '... name_fil: ',trim(adjustl(name_fil(i_filter)))
          write (W3aux,'(e27.8)') PhotFlux(i_filter)
          write (*,'(4x,a,a)') '... Photometric Fluxes: '//trim(adjustl(W3aux))
       end if

       do j_filter=1,Ncalibra
          if ( flux_T01 > 0.0_RP .AND. Mcalibra(j_filter) > -999.0_RP ) then
             MAG_spec(i_filter,j_filter) = mag_Sgal + Mcalibra(j_filter)
          else
             MAG_spec(i_filter,j_filter) = -999.0_RP
          end if
          
          if ( IsShowOn == 1_IB ) then
             write (W1aux,'(e27.8)') MAG_spec(i_filter,j_filter)
             write (W2aux,'(e27.8)') Mcalibra(j_filter)
             
             if ( MAG_spec(i_filter,j_filter) >= 0.0_RP ) then
                W1aux = '+'//trim(adjustl(W1aux))
             end if
             
             if ( Mcalibra(j_filter) >= 0.0_RP ) then
                W2aux = '+'//trim(adjustl(W2aux))
             end if
             
             write (*,'(4x,a,a,a,a)')                                    &
                  '... '//trim(adjustl(Fnamecal(j_filter)))//': ',       &
                  trim(adjustl(W1aux)),                                  &
                  ' Calibration: '//trim(adjustl(W2aux))
          end if
       end do
       
       if ( IsShowOn == 1_IB ) then
          write (*,'(4x,a,a)') '++++++++++++++++++++++++++++++++++++++++++'
       end if
! *** Evaluate magnitudes ***************************************************
    end do

    deallocate( Fnamecal )
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)')   '[EvalFilters]'
    end if

    return
END SUBROUTINE EvalFilters
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_EvalFilters( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_EvalFilters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Sat Oct 27 08:58:18 WEST 2012 ++++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 002                                                          !
!
! 1) EvalFilters
! 2) author_EvalFilters

! ###########################################################################
!     RESUME : Read Filters and retrieve ZERO point calibrations. Original  !
!              routine dates back to 2003-2004.                             !
!                                                                           !
!              Original routines first released:                            !
!                                                                           !
!              Wed Sep 15 09:32:47 WEST 2004                                !
!                                                                           !
!     INPUT    : 01) T_lambda  -> Wavelength                                !
!                02) T_fluxes  -> Transmission curve                        !
!                03) file_dir  -> Filters directory                         !
!                04) file_out  -> In case of error                          !
!                05) Int_Type  -> Integral type                             !
!                06) IsKeepOn  -> Variable:       0 => Run aborted          !
!                07) verbosity -> Optional variable to print & check        !
!                                                                           !
!     OUTPUT   : 01) T_l_Area  -> Area of filter (lambda)                   !
!                02) T_n_Area  -> Area of filter (frequency)                !
!                03) magABsys  -> Calibration for AB system                 !
!                04) magTGsys  -> Calibration for TG system                 !
!                05) standard  -> standard flux for VEGA                    !
!                06) lamb_eff  -> Effective wavelength                      !
!                07) widtheff  -> Effective width                           ! 
!                08) Numb_lbd  -> Number of points in filter                !
!                09) name_fil  -> Name of bandpass/filter                   !
!                10) Nfilters  -> Number of filters                         !
!                                                                           !
!     EXTRA ROUTINES : EvalTransmission and IntegralALL                     !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!     Checked: Mon Sep 25 04:06:38 PM WEST 2023                             !
!                                                                           !
!     LOG: 25/09/2023: Correction for flux_T01 > 0.0 in magABsys            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PropFilters( T_lambda,T_fluxes,Ntlambda,Nfilters,T_l_Area,       &
                        T_n_Area,magABsys,magTGsys,standard,lamb_eff,       &
                        widtheff,Numb_lbd,lambVega,fluxVega,NVegalbd,       &
                        lamb_Sun,flux_Sun,NSun_lbd,lamb_FBD,flux_FBD,       &
                        NFBD_lbd,IsKeepOn,Int_Type,verbosity )
    use ModDataType
    implicit none
    real     (kind=RP), parameter :: cl_speed=299792.4580_RP
    real     (kind=RP), parameter :: fac4Pid2=1.19649518e40_RP
    
    integer  (kind=IB), intent(in) :: Ntlambda,Nfilters
    integer  (kind=IB), intent(in) :: NVegalbd,NSun_lbd,NFBD_lbd
    integer  (kind=IB), intent(out) :: IsKeepOn

    integer  (kind=IB), optional :: Int_Type,verbosity
    integer  (kind=IB) :: IsShowOn,i_filter,N_lambda,IntegraT,IsOFFNum
    integer  (kind=IB), dimension(Nfilters), intent(in) :: Numb_lbd

    real     (kind=RP), target, dimension(Ntlambda,Nfilters),intent(in) ::  &
                                                                  T_lambda, &
                                                                  T_fluxes
    real     (kind=RP), pointer, dimension(:) :: P1lambda,P1fluxes
    real     (kind=RP), dimension(Nfilters), intent(out) :: T_l_Area,       &
                                                            T_n_Area,       &
                                                            magABsys,       &
                                                            magTGsys,       &
                                                            standard,       &
                                                            lamb_eff,       &
                                                            widtheff

    real     (kind=RP), dimension(NVegalbd), intent(in) :: lambVega,fluxVega
    real     (kind=RP), dimension(NSun_lbd), intent(in) :: lamb_Sun,flux_Sun
    real     (kind=RP), dimension(NFBD_lbd), intent(in) :: lamb_FBD,flux_FBD
    real     (kind=RP), dimension(:), allocatable :: aux_flux

    real     (kind=RP) :: lambda_i,lambda_f,aux_area,mAB_Vega,mTG_Vega,     &
                          flux_T01,flux_T02,flux_T03

    character (len=CH) :: W1aux

    !f2py real     (kind=RP), intent(in)  :: T_lambda
    !f2py real     (kind=RP), intent(in)  :: T_fluxes
    !f2py real     (kind=RP), intent(out) :: T_l_Area
    !f2py real     (kind=RP), intent(out) :: T_n_Area
    !f2py real     (kind=RP), intent(out) :: magABsys
    !f2py real     (kind=RP), intent(out) :: magTGsys
    !f2py real     (kind=RP), intent(out) :: standard
    !f2py real     (kind=RP), intent(out) :: lamb_eff
    !f2py real     (kind=RP), intent(out) :: widtheff
    !f2py integer  (kind=IB), intent(in)  :: Numb_lbd
    !f2py integer  (kind=IB), intent(hide), depend(T_lambda) :: Ntlambda=shape(T_lambda,0), Nfilters=shape(T_lambda,1)

    !f2py real     (kind=RP), intent(in) :: lambVega
    !f2py real     (kind=RP), intent(in) :: fluxVega
    !f2py integer  (kind=IB), intent(hide), depend(lambVega) :: NVegalbd=shape(lambVega,0)  
    
    !f2py real     (kind=RP), intent(in) :: lamb_Sun
    !f2py real     (kind=RP), intent(in) :: flux_Sun
    !f2py integer  (kind=IB), intent(hide), depend(lamb_Sun) :: NSun_lbd=shape(lamb_Sun,0)  
    
    !f2py real     (kind=RP), intent(in) :: lamb_FBD
    !f2py real     (kind=RP), intent(in) :: flux_FBD
    !f2py integer  (kind=IB), intent(hide), depend(lamb_FBD) :: NFBD_lbd=shape(lamb_FBD,0)  
    
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    
    !f2py integer  (kind=IB), optional :: Int_Type=2
    !f2py integer  (kind=IB), optional :: verbosity=0
    
    interface
       subroutine EvalTransmission( L_lambda,S_fluxes,N_lambda,T_lambda,   &
                                    T_fluxes,Ntlambda,Int_Type,fluxtran,   &
                                    IsKeepOn,verbosity )
          use ModDataType
          implicit none
          integer  (kind=IB), intent(out) :: IsKeepOn
          integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
          integer  (kind=IB), optional :: verbosity
          real     (kind=RP), dimension(N_lambda), intent(in) :: L_lambda,  &
                                                                 S_fluxes
          real     (kind=RP), dimension(Ntlambda), intent(in) :: T_lambda,  &
                                                                 T_fluxes
          real     (kind=RP), intent(out) :: fluxtran
        end subroutine EvalTransmission
        subroutine EvalTFluxes( O_lambda,O_fluxes,Nspeclbd,T_lambda,        &
                                T_fluxes,N_lambda,T_l_area,fluxtran,        &
                                lamb_eff,IskeepOn,verbosity )
            use ModDataType
            integer  (kind=IB), intent(out) :: IsKeepOn
            integer  (kind=IB), intent(in) :: N_lambda,Nspeclbd
            integer  (kind=IB), optional :: verbosity
            real     (kind=RP), dimension(Nspeclbd), intent(in) :: O_lambda,&
                                                                   O_fluxes
            real     (kind=RP), dimension(N_lambda), intent(in) :: T_lambda,&
                                                                   T_fluxes
            real     (kind=RP), intent(out) :: lamb_eff,fluxtran
            real     (kind=RP), intent(in) :: T_l_area
        end subroutine EvalTFluxes
        real (kind=RP) function IntegralALL( SXvalues,SYvalues,lambda_i,    &
                                             lambda_f,N_lambda,IsKeepOn,    &
                                             Int_Type,verbosity )
            use ModDataType
            integer  (kind=IB), intent(out) :: IsKeepOn
            integer  (kind=IB), intent(in) :: N_lambda,Int_Type
            integer  (kind=IB), optional :: verbosity
            real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,&
                                                                   SYvalues
            real     (kind=RP), intent(in) :: lambda_i,lambda_f
         end function IntegralALL
         real (kind=RP) function lamb_effective( T_lambda,T_fluxes,Ntlambda,&
                                                 L_lambda,S_fluxes,N_lambda,&
                                                 IskeepOn,Int_Type )
           use ModDataType
           implicit none
           integer  (kind=IB), intent(out) :: IsKeepOn
           integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
           real     (kind=RP), intent(in), dimension(Ntlambda) :: T_lambda, &
                                                                  T_fluxes
           real     (kind=RP), intent(in), dimension(N_lambda) :: L_lambda, &
                                                                  S_fluxes
         end function lamb_effective
    end interface

    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    if ( present(Int_Type) ) then
       IntegraT = Int_Type
    else
       IntegraT= 2_IB
    end if
    
!  *** Reset Filter-STUFF ***************************************************
    T_l_Area = -999.0_RP
    T_n_Area = -999.0_RP
    magABsys = -999.0_RP
    magTGsys = -999.0_RP
    standard = -999.0_RP
    lamb_eff = -999.0_RP

    flux_T01 = -999.0_RP
    flux_T02 = -999.0_RP
    flux_T03 = -999.0_RP

    IsOFFNum = +00000_IB
!  *** Reset Filter-STUFF ***************************************************

    !ilastnum = index(file_dir,' ') - 1

    allocate ( aux_flux(Ntlambda) )

    ! Ancient part of the code, where it read the content
    
! *** Read calibration stars ************************************************
!     RESUME : VEGA spectrum.                                               !
!              Intrinsic Flux: erg/s/cm^2/A                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !arq_fil1 = file_dir(1:ilastnum)//'VegaLR.dat'
    !open  (21,status='old',file=arq_fil1,ERR=22)
    !read  (21,*,ERR=22) arq_lixo,NVegalbd
    !
    !allocate ( lambVega(NVegalbd) )
    !allocate ( fluxVega(NVegalbd) )
    !
    !do i_lambda=1,NVegalbd
    !    read  (21,*,ERR=22) lambVega(i_lambda),fluxVega(i_lambda)
    !end do
    !close (21)

! *** SUN spectrum **********************************************************
!     Intrinsic flux: erg/s/A.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !arq_fil1 = file_dir(1:ilastnum)//'Sun_LR.dat'
    !open(21,status='old',file=arq_fil1,ERR=22)
    !read(21,*,ERR=22) arq_lixo,NSun_lbd
    !
    !allocate ( lamb_Sun(NSun_lbd) )
    !allocate ( flux_Sun(NSun_lbd) )
    !
    !do i_lambda=1,NSun_lbd
    !    read(21,*,ERR=22) lamb_Sun(i_lambda),flux_Sun(i_lambda)
    !end do
    !close(21)

! *** Reading of the spectrum of BD+17o4708 *********************************
!     RESUME : F subdwarf used to calibrate the Thuan and Gunn system.      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !arq_fil1 = file_dir(1:ilastnum)//'BD+17o4708.dat'
    !arq_fil1 = file_dir(1:ilastnum)//'BD+17d4708.dat'
    !open(21,status='old',file=arq_fil1,ERR=22)
    !read(21,*,ERR=22) arq_lixo,NFBD_lbd
    !
    !allocate ( lamb_FBD(NFBD_lbd) )
    !allocate ( flux_FBD(NFBD_lbd) )
    !
    !do i_lambda=1,NFBD_lbd
    !    read(21,*,ERR=22) lamb_FBD(i_lambda),flux_FBD(i_lambda)
    !end do
    !close(21)
! *** Read calibration stars ************************************************

! *** Compute in Filters ****************************************************
    do i_filter=1,Nfilters
       T_l_Area(i_filter) = 0.0_RP
       T_n_Area(i_filter) = 0.0_RP

       N_lambda = Numb_lbd(i_filter)
       lambda_i = T_lambda(00000001,i_filter)
       lambda_f = T_lambda(N_lambda,i_filter)

       if ( N_lambda == 1_IB ) then
          !# If this happens it should be 1.0, but be careful!
          T_l_Area(i_filter) = 1.0_RP !* T_fluxes(1,i_filter)
          T_n_Area(i_filter) = 1.0_RP !* T_fluxes(1,i_filter)
          
          !flux_T01 = ?
          !flux_T02 = ?
          !flux_T03 = ?
          widtheff(i_filter) = 1.000_RP
          lamb_eff(i_filter) = lambda_i
          standard(i_filter) = flux_T01
       else
          if (associated (P1lambda)) nullify( P1lambda )
          if (associated (P1fluxes)) nullify( P1fluxes )
          
          P1lambda => T_lambda(1:N_lambda,i_filter)
          P1fluxes => T_fluxes(1:N_lambda,i_filter)

          T_l_Area(i_filter) = IntegralALL( P1lambda,P1fluxes,lambda_i,     &
                                            lambda_f,N_lambda,IsKeepOn,     &
                                            IntegraT,IsOFFNum )

          aux_area = T_l_Area(i_filter)
          
          aux_flux(1:N_lambda) = P1fluxes(1:N_lambda)                       &
                               / (P1lambda(1:N_lambda)                      &
                               * P1lambda(1:N_lambda)) * cl_speed
          aux_flux(1:N_lambda) = aux_flux(1:N_lambda) !* 1.0e13_RP
          
          T_n_Area(i_filter) = IntegralALL( P1lambda,aux_flux,lambda_i,     &
                                            lambda_f,N_lambda,IsKeepOn,     &
                                            IntegraT,IsOFFNum )
          T_n_Area(i_filter) = T_n_Area(i_filter) * 1.0e13_RP

! *** Evaluate transmission flux ********************************************
          !call EvalTFluxes( lambVega,fluxVega,P1lambda,P1fluxes,aux_area,  &
          !                  N_lambda,flux_T01,lamb_T01,NVegalbd,IskeepOn,  &
          !                  IsShowOn )
          !call EvalTFluxes( lamb_Sun,flux_Sun,P1lambda,P1fluxes,aux_area,  &
          !                  N_lambda,flux_T02,lamb_T02,NSun_lbd,IskeepOn,  &
          !                  IsShowOn )
          !call EvalTFluxes( lamb_FBD,flux_FBD,P1lambda,P1fluxes,aux_area,  &
          !                  N_lambda,flux_T03,lamb_T03,NFBD_lbd,IskeepOn,  &
          !                  IsShowOn )

          call EvalTransmission( lambVega,fluxVega,NVegalbd,P1lambda,       &
                                 P1fluxes,N_lambda,IntegraT,flux_T01,       &
                                 IsKeepOn,IsOFFNum )
          call EvalTransmission( lamb_Sun,flux_Sun,NSun_lbd,P1lambda,       &
                                 P1fluxes,N_lambda,IntegraT,flux_T02,       &
                                 IsKeepOn,IsOFFNum )
          call EvalTransmission( lamb_FBD,flux_FBD,NFBD_lbd,P1lambda,       &
                                 P1fluxes,N_lambda,IntegraT,flux_T03,       &
                                 IsKeepOn,IsOFFNum )

          !write (*,*) standard(i_filter),flux_T01
          standard(i_filter) = flux_T01
          lamb_eff(i_filter) = lamb_effective( P1lambda,P1fluxes,N_lambda,  &
                                               lambVega,fluxVega,NVegalbd,  &
                                               IskeepOn,IntegraT )

          ! Effective width
          ! width = integral( T x dl) / max(T)
          if ( maxval( P1fluxes ) > 0.0_RP ) then
             widtheff(i_filter) = T_l_Area(i_filter) / maxval( P1fluxes )
          else
             widtheff(i_filter) = -999.0_RP
          end if
! *** Evaluate transmission flux ********************************************
       end if

! *** Magnitudes AB system **************************************************
       if ( flux_T01 > 0.0_RP ) then
          mAB_Vega = -2.5_RP                                                &
                * log10( flux_T01*T_l_area(i_filter)/T_n_area(i_filter) )   &
                - 48.60_RP       
       else
          mAB_Vega = -999.0_RP
       end if
       magABsys(i_filter) = mAB_Vega       
! *** Magnitudes AB system **************************************************

! *** Magnitudes Thuan and Gunn system **************************************
       if ( flux_T03 > 0.0_RP ) then
          mTG_Vega = -2.50_RP * log10(flux_T01/flux_T03) + 9.5_RP
       else
          mTG_Vega = -999.0_RP
       end if
       magTGsys(i_filter) = mTG_Vega
! *** Magnitudes Thuan and Gunn system **************************************

! *** Print screen **********************************************************       
       if ( IsShowOn == 1_IB ) then
          write (W1aux,'(f18.12)') T_l_Area(i_filter)
          write (*,'(4x,a,a)') '... T_l_Area: ',trim(adjustl(W1aux))
          write (W1aux,'(e18.12)') T_n_Area(i_filter)
          write (*,'(4x,a,a)') '... T_n_Area: ',trim(adjustl(W1aux))
          write (W1aux,'(f18.12)') mAB_Vega
          write (*,'(4x,a,a)') '... M_ABVega: ',trim(adjustl(W1aux))
          write (W1aux,'(f18.12)') mTG_Vega
          write (*,'(4x,a,a)') '... T_G_Vega: ',trim(adjustl(W1aux))
          write (W1aux,'(f18.12)') lamb_eff(i_filter)
          write (*,'(4x,a,a)') '... lamb_eff: ',trim(adjustl(W1aux))
          write (W1aux,'(f18.12)') widtheff(i_filter)
          write (*,'(4x,a,a)') '... widtheff: ',trim(adjustl(W1aux))
          write (*,'(4x,a,a)') '++++++++++++++++++++++++++++++++++++++++++'
       end if
        
    end do
! *** Compute in Filters ****************************************************

    deallocate ( aux_flux )

    if ( IsShowOn == 1_IB ) then
        write (*,'(4x,a)')   '[PropFilters]'
    end if

    return

! *** Issue warning - NO READING OF SPECTRUM ********************************
24  format(4x,65('?'))
!22  write (*,24)
    write (*,24)
    write (*,'(4x,a)') '[ReadFilters] WARNING! NO SPECTRUM READ @@@@@@@'
    write (*,24)

!23  write (*,24)
    write (*,24)
    write (*,'(4x,a)') '[ReadFilters] WARNING! PROBLEM - Filters List @'
    write (*,24)

    IsKeepOn = 0_IB
! *** Issue warning - NO EXTRAPOLATION **************************************

END SUBROUTINE PropFilters
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_PropFilters( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_PropFilters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Compute effective lambda and flux convoluted with the        !
!              transmission curve.                                          !
!                                                                           !
!     INPUT    : 01) O_lambda -> Wavelength                                 !
!                02) O_fluxes -> Fluxes                                     !
!                03) T_lambda -> Transmission lambda                        !
!                04) T_fluxes -> Transmission                               !
!                05) T_l_area -> Area computed with lambda                  !
!                06) N_lambda -> Number of filter points                    !
!                07) Nspeclbd -> Number of spectrum points                  !
!                08) IsKeepOn -> Variable:       0 => Run aborted           !
!                09) verbosity -> Optional variable to print and check      !
!                                                                           !
!     OUTPUT   : 01) fluxtran -> Flux inside the filter                     !
!                02) lamb_eff -> Effective lambda                           !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvalTFluxes( O_lambda,O_fluxes,Nspeclbd,T_lambda,T_fluxes,       &
                        N_lambda,T_l_area,fluxtran,lamb_eff,IskeepOn,       &
                        verbosity )
    use ModDataType
    implicit none
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: N_lambda,Nspeclbd
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,jnum_aux,knum_aux
    real     (kind=RP), dimension(Nspeclbd), intent(in) :: O_lambda,O_fluxes
    real     (kind=RP), dimension(N_lambda), intent(in) :: T_lambda,T_fluxes
    real     (kind=RP), intent(out) :: lamb_eff,fluxtran
    real     (kind=RP), intent(in) :: T_l_area
    real     (kind=RP) :: lambda_i,lambda_f,transinf,flux_inf,transsup,     &
                          flux_sup

    !f2py real     (kind=RP), intent(in)  :: O_lambda
    !f2py real     (kind=RP), intent(in)  :: O_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(O_lambda) :: Nspeclbd=shape(O_lambda,0)

    !f2py real     (kind=RP), intent(in)  :: T_lambda
    !f2py real     (kind=RP), intent(in)  :: T_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(T_lambda) :: Ntlambda=shape(T_lambda,0)

    !f2py real     (kind=RP), intent(in)  :: T_l_Area
    !f2py real     (kind=RP), intent(out) :: fluxtran
    !f2py real     (kind=RP), intent(out) :: lamb_eff
    
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    
    !f2py integer  (kind=IB), optional :: verbosity=0
    
    interface
        real (kind=RP) function SLin_interp( xorg_ini,xorg_fin,yorg_ini,    &
                                             yorg_fin,newvalue )
            use ModDataType
            real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,   &
                                              yorg_fin, newvalue
        end function SLin_interp
        real (kind=RP) function SLog_interp( xorg_ini,xorg_fin,yorg_ini,    &
                                             yorg_fin,newvalue )
            use ModDataType
            real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,   &
                                              yorg_fin,newvalue
        end function SLog_interp
    end interface

    IsKeepOn = 1_IB

    if ( verbosity == 1_IB ) then
        IsShowOn = 1_IB
    end if

    lamb_eff = 0.0_RP
    fluxtran = 0.0_RP

    if ( O_lambda(1) > T_lambda(1) .OR.                                     &
         O_lambda(Nspeclbd) < T_lambda(N_lambda) ) then
        fluxtran = -999.000_RP
    else
        knum_aux = 1
        do while( O_lambda(knum_aux+1) < T_lambda(1) )
            knum_aux = knum_aux + 1
        end do
        jnum_aux = 1
        lambda_i = T_lambda(jnum_aux)
        lambda_f = lambda_i
        transinf = 0.0_RP
        flux_inf = 0.0_RP

        do while ( jnum_aux+1 <= N_lambda .AND. knum_aux+1 <= Nspeclbd )
            if ( T_lambda(jnum_aux+1) < O_lambda(knum_aux+1) ) then
                lambda_f = T_lambda(jnum_aux+1)
                flux_sup = SLog_interp( O_lambda(knum_aux),                 &
                                        O_lambda(knum_aux+1),               &
                                        O_fluxes(knum_aux),                 &
                                        O_fluxes(knum_aux+1),lambda_f )
                transsup = T_fluxes(jnum_aux+1)! * lambda_f**typetrans(i)
                jnum_aux = jnum_aux+1
            else
                lambda_f = O_lambda(knum_aux+1)
                flux_sup = O_fluxes(knum_aux+1)
                transsup = SLin_interp( T_lambda(jnum_aux),                 &
                                        T_lambda(jnum_aux+1),               &
                                        T_fluxes(jnum_aux),                 &
                                        T_fluxes(jnum_aux+1),               &
                                        lambda_f )
                transsup = transsup !* lambda_f**typetrans(i)
                knum_aux = knum_aux+1
            end if
             
            fluxtran = fluxtran + (transinf*flux_inf+transsup*flux_sup)     &
                                * (lambda_f-lambda_i) / 2.0_RP
            lamb_eff = lamb_eff + (transinf*flux_inf+transsup*flux_sup)     &
                                * (lambda_f-lambda_i)/2.0_RP                &
                     * (lambda_f+lambda_i)/2.0_RP
            lambda_i = lambda_f
            flux_inf = flux_sup
            transinf = transsup
        end do
        lamb_eff = lamb_eff / fluxtran
        fluxtran = fluxtran / T_l_area
     end if

     return
END SUBROUTINE EvalTFluxes
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_EvalTFluxes( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_EvalTFluxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Compute flux convoluted with the transmission curve.         !
!                                                                           !
!     INPUT    : 01) O_lambda  -> Wavelength                                !
!                02) O_fluxes  -> Fluxes                                    !
!                03) T_lambda  -> Transmission lambda                       !
!                04) T_fluxes  -> Transmission                              !
!                05) T_l_area  -> Area computed with lambda                 !
!                06) N_lambda  -> Number of filter points                   !
!                07) Nspeclbd  -> Number of spectrum points                 !
!                08) IsKeepOn  -> Variable:       0 => Run aborted          !
!                09) verbosity -> Optional variable to print and check      !
!                                                                           !
!     OUTPUT   : 01) fluxtran  -> Flux inside the filter                    !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvalTransmission( L_lambda,S_fluxes,N_lambda,T_lambda,T_fluxes,  &
                             Ntlambda,Int_Type,fluxtran,IsKeepOn,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,Nlbd_aux,ilastval,IsOFFNum,ind1,ind2

    real     (kind=RP), dimension(N_lambda), intent(in) :: L_lambda,        &
                                                           S_fluxes
    real     (kind=RP), dimension(Ntlambda), intent(in) :: T_lambda,        &
                                                           T_fluxes
    real     (kind=RP), allocatable, dimension(:) :: aux__lbd,aux_flux,     &
                                                     aux__fil
    real     (kind=RP), intent(out) :: fluxtran
    real     (kind=RP) :: lambda_i,lambda_f,dtlambda,d_lambda

    !f2py real     (kind=RP), intent(in)  :: L_lambda
    !f2py real     (kind=RP), intent(in)  :: S_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(L_lambda) :: N_lambda=shape(L_lambda,0)

    !f2py real     (kind=RP), intent(in)  :: T_lambda
    !f2py real     (kind=RP), intent(in)  :: T_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(T_lambda) :: Ntlambda=shape(T_lambda,0)

    !f2py real     (kind=RP), intent(out) :: fluxtran

    !f2py integer  (kind=IB), intent(in)  :: In_Type
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    
    !f2py  integer  (kind=IB), optional :: verbosity=0
    
    interface
       subroutine LINinterpol( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,&
                               nold_vec,ilastval,IsKeepOn,verbosity )

         use ModDataType
         integer  (kind=IB), intent(in out) :: ilastval
         integer  (kind=IB), intent(in) :: nold_vec, nxyvalue
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(0:nold_vec-1), intent(in) ::         &
                                                                   xold_vec,&
                                                                   yold_vec
         real     (kind=RP), dimension(0:nxyvalue-1), intent(in) :: xx_value
         real     (kind=RP), dimension(0:nxyvalue-1), intent(out) :: yy_value
       end subroutine LINinterpol
    
       real (kind=RP) function IntegralALL( SXvalues,SYvalues,lambda_i,     &
                                            lambda_f,N_lambda,IsKeepOn,     &
                                            Int_Type,verbosity )
         use ModDataType
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: N_lambda,Int_Type
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,   &
                                                                SYvalues
         real     (kind=RP), intent(in) :: lambda_i,lambda_f
       end function IntegralALL
    end interface

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if

    IsOFFNum = 0_IB

    ! Interpolate to lowest dlambda
    dtlambda = minval( array=T_lambda(2:Ntlambda)-T_lambda(1:Ntlambda-1) )

    if ( L_lambda(1) <= T_lambda(1) .AND. L_lambda(N_lambda) >= T_lambda(Ntlambda) ) then
       d_lambda = minval( array=L_lambda(2:N_lambda)-L_lambda(1:N_lambda-1),   &
            mask=L_lambda(2:N_lambda) >= T_lambda(1)                .AND.      &
                 L_lambda(2:N_lambda) <= T_lambda(Ntlambda) )
    else
       d_lambda = -999.0_RP
       fluxtran = -999.0_RP
       IsKeepOn = 0_IB
       return
    end if
    
    !write (*,*) dtlambda,d_lambda
    
    ilastval = -999_IB
    if ( dtlambda < d_lambda ) then
       Nlbd_aux = Ntlambda
       allocate( aux__lbd(Ntlambda) )
       allocate( aux_flux(Ntlambda) )
       
       aux__lbd(1:Ntlambda) = T_lambda(1:Ntlambda)
       call LINinterpol( aux__lbd,aux_flux,Nlbd_aux,L_lambda,S_fluxes,N_lambda,ilastval,IsKeepOn,IsOFFNum )
       
       lambda_i = T_lambda(1)
       lambda_f = T_lambda(Ntlambda)

       !flux = integrate[(f_lambda*lambda*T(lambda*(1+z))*dlambda)]
       !norm = integrate[(lambda*T(lambda*(1+z))*dlambda))]
       !photometry = flux/norm
        
       fluxtran = IntegralALL( aux__lbd,aux__lbd*T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
       fluxtran = fluxtran / IntegralALL( aux__lbd,aux__lbd*T_fluxes,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
       
       deallocate( aux__lbd )
       deallocate( aux_flux )
    else
       ind1 = maxloc( array=L_lambda, mask=L_lambda <= T_lambda(1), dim=1 )
       ind2 = minloc( array=L_lambda, mask=L_lambda >= T_lambda(Ntlambda), dim=1 )

       !write (*,*) ind1,ind2,L_lambda(ind1),L_lambda(ind2),S_fluxes(ind1),S_fluxes(ind2)

       if ( ind1 > 0 .AND. ind2 > 0 .AND. ind1 < N_lambda .AND. ind2 < N_lambda .AND. ind1 <= ind2 ) then
          aux__lbd = pack( array=L_lambda, mask=L_lambda >= L_lambda(ind1) .AND. L_lambda <= L_lambda(ind2)  )
          aux_flux = pack( array=S_fluxes, mask=L_lambda >= L_lambda(ind1) .AND. L_lambda <= L_lambda(ind2)  )
          Nlbd_aux = size( aux__lbd )

          !write(*,*) size(aux__lbd),size(aux_flux),Nlbd_aux
          
          !write (*,*) Nlbd_aux,T_lambda(1),T_lambda(Ntlambda),aux__lbd(1),aux__lbd(Nlbd_aux)
          if ( Nlbd_aux > 1_IB  ) then
             allocate( aux__fil(Nlbd_aux) )

             lambda_i = T_lambda(1) !aux__lbd(1)
             lambda_f = T_lambda(Ntlambda) !aux__lbd(Nlbd_aux)
             
             call LINinterpol( aux__lbd,aux__fil,Nlbd_aux,T_lambda,T_fluxes,Ntlambda,ilastval,IsKeepOn,IsOFFNum )
             where ( aux__fil(1:Nlbd_aux) < 0.0_RP )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where
             where ( aux__lbd(1:Nlbd_aux) < lambda_i )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where
             where ( aux__lbd(1:Nlbd_aux) > lambda_f )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where

             !write (*,*) size(aux__lbd),Nlbd_aux,size(aux__fil),size(aux_flux)
             
             fluxtran = IntegralALL( aux__lbd,aux__lbd*aux_flux*aux__fil,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
             fluxtran = fluxtran / IntegralALL( aux__lbd,aux__lbd*aux__fil,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )

             !write (*,*) fluxtran
             deallocate( aux__fil )
          end if
          deallocate( aux__lbd )
          deallocate( aux_flux )
       end if
    end if
    
    return
END SUBROUTINE EvalTransmission
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_EvalTransmission( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_EvalTransmission
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Linear interpolation given 2 points (x1,y1) and (x2,y2).     !
!                                                                           !
!     INPUT    : 01) xorg_ini -> x1                                         !
!                02) xorg_fin -> x2                                         !
!                03) yorg_ini -> y1                                         !
!                04) yorg_fin -> y2                                         !
!                05) newvalue -> new x value                                !
!                                                                           !
!     OUTPUT   : 01) SLin_interp -> new y value                             !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SLin_interp( xorg_ini,xorg_fin,yorg_ini,yorg_fin,   &
                                     newvalue )
    use ModDataType
    implicit none
    real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,yorg_fin,  &
                                      newvalue

    !f2py real     (kind=RP), intent(in)   :: xorg_ini
    !f2py real     (kind=RP), intent(in)   :: xorg_fin
    !f2py real     (kind=RP), intent(in)   :: yorg_ini
    !f2py real     (kind=RP), intent(in)   :: yorg_fin
    !f2py real     (kind=RP), intent(in)   :: newvalue
    !f2py real     (kind=RP), intent(out)  :: SLin_interp
    
    SLin_interp = yorg_ini + (yorg_fin-yorg_ini) * (newvalue-xorg_ini)      &
                / (xorg_fin-xorg_ini)
    return
END FUNCTION SLin_interp
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_SLin_interp( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_SLin_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Logarithmically interpolate given 2 points (x1,y1) and       !
!              (x2,y2).                                                     !
!                                                                           !
!     INPUT    : 01) xorg_ini -> x1                                         !
!                02) xorg_fin -> x2                                         !
!                03) yorg_ini -> y1                                         !
!                04) yorg_fin -> y2                                         !
!                05) newvalue -> new x value                                !
!                                                                           !
!     OUTPUT   : 01) SLog_interp -> new y value                             !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SLog_interp( xorg_ini,xorg_fin,yorg_ini,yorg_fin,   &
                                     newvalue )
    use ModDataType
    implicit none
    real     (kind=RP), parameter  :: setlimit=1.0e-37_RP
    real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,yorg_fin,  &
                                      newvalue

    !f2py real     (kind=RP), intent(in)   :: xorg_ini
    !f2py real     (kind=RP), intent(in)   :: xorg_fin
    !f2py real     (kind=RP), intent(in)   :: yorg_ini
    !f2py real     (kind=RP), intent(in)   :: yorg_fin
    !f2py real     (kind=RP), intent(in)   :: newvalue
    !f2py real     (kind=RP), intent(out)  :: SLog_interp
    
    if ( yorg_ini > setlimit .AND. yorg_fin > setlimit ) then
       SLog_interp = yorg_ini                                               &
                   * (yorg_fin/yorg_ini)**(log(newvalue/xorg_ini)           &
                   / log(xorg_fin/xorg_ini))
    else
        SLog_interp = 0.0_RP
    endif

    return
END FUNCTION SLog_interp
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_SLog_interp( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_SLog_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Effective wavelength from filter:                            !
!                                                                           !
!              lamb_effective =                                             !
!                                                                           !
!              integral(lamb * T * Vega dlamb) / integral(T * Vega dlamb)   !
!                                                                           !
!     INPUT    : 01) T_lambda -> Wavelength from transmission filter        !
!                02) T_fluxes -> Transmission values                        !
!                03) Ntlambda -> # of points in transmission filter         !
!                04) L_lambda -> Wavelength of calibration star (e.g. Vega) !
!                05) S_fluxes -> "Flux" of calibration star (e.g. Vega)     !
!                06)                                                        !
!                                                                           !
!     OUTPUT   : 01) SLog_interp -> new y value                             !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : Comparison of legacy fortran formula with the one from       !
!              pyphot, for the release of the code. But it is a more        !
!              general function.                                            !
!                                                                           !
!     Checked: 2021                                                         !
!                                                                           !
!   """Unitwise Effective wavelength                                        !
!   leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)                !
!   """                                                                     !
!   with Vega() as v:                                                       !
!       s = self.reinterp(v.wavelength)                                     !
!       w = s._wavelength                                                   !
!       leff = np.trapz(w * s.transmit * v.flux.magnitude, w, axis=-1)      !
!       leff /= np.trapz(s.transmit * v.flux.magnitude, w, axis=-1)         !
!   if self.wavelength_unit is not None:                                    !
!       return leff * unit[self.wavelength_unit]                            !
!   else:                                                                   !
!       return leff                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION lamb_effective( T_lambda,T_fluxes,Ntlambda,L_lambda,&
                                        S_fluxes,N_lambda,IsKeepOn,Int_Type )
    use ModDataType
    implicit none
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
    integer  (kind=IB) :: Nlbd_aux,ilastval,IsOFFNum,ind1,ind2
    
    real     (kind=RP), intent(in), dimension(Ntlambda) :: T_lambda,T_fluxes
    real     (kind=RP), intent(in), dimension(N_lambda) :: L_lambda,S_fluxes

    real     (kind=RP), allocatable, dimension(:) :: aux__lbd,aux_flux,aux__fil

    real     (kind=RP) :: dtlambda,d_lambda,lambda_i,lambda_f

    !f2py real     (kind=RP), intent(in)  :: T_lambda
    !f2py real     (kind=RP), intent(in)  :: T_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(T_lambda) :: Ntlambda=shape(T_lambda,0)

    
    !f2py real     (kind=RP), intent(in)  :: L_lambda
    !f2py real     (kind=RP), intent(in)  :: S_fluxes
    !f2py integer  (kind=IB), intent(hide), depend(L_lambda) :: N_lambda=shape(L_lambda,0)

    !f2py real     (kind=RP), intent(out) :: lamb_effective

    !f2py integer  (kind=IB), intent(in)  :: In_Type
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    
    interface
       subroutine LINinterpol( xx_value,yy_value,nxyvalue,xold_vec,yold_vec, &
                               nold_vec,ilastval,IsKeepOn,verbosity )

         use ModDataType

         implicit none
         integer  (kind=IB), intent(in out) :: ilastval
         integer  (kind=IB), intent(in) :: nold_vec, nxyvalue
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(0:nold_vec-1), intent(in) :: xold_vec,&
                                                                    yold_vec
         real     (kind=RP), dimension(0:nxyvalue-1), intent(in) :: xx_value
         real     (kind=RP), dimension(0:nxyvalue-1), intent(out) :: yy_value
       end subroutine LINinterpol
    
       real (kind=RP) function IntegralALL( SXvalues,SYvalues,lambda_i,      &
                                            lambda_f,N_lambda,IsKeepOn,      &
                                            Int_Type,verbosity )
         use ModDataType
         integer  (kind=IB), intent(in out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: N_lambda,Int_Type
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,    &
                                                                SYvalues
         real     (kind=RP), intent(in) :: lambda_i,lambda_f
       end function IntegralALL
    end interface

    IsOFFNum = 0_IB
    lamb_effective = 0.0_RP
    
    ! Interpolate to lowest dlambda
    dtlambda = minval( T_lambda(2:Ntlambda)-T_lambda(1:Ntlambda-1) )
    d_lambda = minval( L_lambda(2:N_lambda)-L_lambda(1:N_lambda-1), &
                       mask=L_lambda(2:N_lambda) >= T_lambda(1) .AND. L_lambda(2:N_lambda) <= T_lambda(Ntlambda)  )

    !write (*,*) dtlambda,d_lambda
    
    ilastval = -999_IB
    if ( dtlambda < d_lambda ) then
       Nlbd_aux = Ntlambda
       allocate( aux__lbd(Ntlambda) )
       allocate( aux_flux(Ntlambda) )
       
       aux__lbd(1:Ntlambda) = T_lambda(1:Ntlambda)
       call LINinterpol( aux__lbd,aux_flux,Nlbd_aux,L_lambda,S_fluxes,N_lambda,ilastval,IsKeepOn,IsOFFNum )
       
       lambda_i = T_lambda(1)
       lambda_f = T_lambda(Ntlambda)
       lamb_effective = IntegralALL( aux__lbd,aux__lbd*T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
       lamb_effective = lamb_effective &
                      / IntegralALL( aux__lbd,T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
       
       deallocate( aux__lbd )
       deallocate( aux_flux )
    else
       ind1 = maxloc( array=L_lambda, mask=L_lambda <= T_lambda(1), dim=1 )
       ind2 = minloc( array=L_lambda, mask=L_lambda >= T_lambda(Ntlambda), dim=1 )

       !write (*,*) ind1,ind2,L_lambda(ind1),L_lambda(ind2)

       if ( ind1 > 0 .AND. ind2 > 0 .AND. ind1 < N_lambda .AND. ind2 < N_lambda .AND. ind1 <= ind2 ) then
          aux__lbd = pack( array=L_lambda, mask=L_lambda >= L_lambda(ind1) .AND. L_lambda <= L_lambda(ind2)  )
          aux_flux = pack( array=S_fluxes, mask=L_lambda >= L_lambda(ind1) .AND. L_lambda <= L_lambda(ind2)  )
          Nlbd_aux = size( aux__lbd )

          !write(*,*) size(aux__lbd),size(aux_flux),Nlbd_aux
          
          !write (*,*) Nlbd_aux,T_lambda(1),T_lambda(Ntlambda),aux__lbd(1),aux__lbd(Nlbd_aux)
          if ( Nlbd_aux > 1_IB  ) then
             allocate( aux__fil(Nlbd_aux) )

             lambda_i = T_lambda(1) !aux__lbd(1)
             lambda_f = T_lambda(Ntlambda) !aux__lbd(Nlbd_aux)
             
             call LINinterpol( aux__lbd,aux__fil,Nlbd_aux,T_lambda,T_fluxes,Ntlambda,ilastval,IsKeepOn,IsOFFNum )
             where ( aux__fil(1:Nlbd_aux) < 0.0_RP )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where
             where ( aux__lbd(1:Nlbd_aux) < lambda_i )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where
             where ( aux__lbd(1:Nlbd_aux) > lambda_f )
                aux__fil(1:Nlbd_aux) = 0.0_RP
             end where

             !write (*,*) size(aux__lbd),Nlbd_aux,size(aux__fil),size(aux_flux)
             
             lamb_effective = IntegralALL( aux__lbd,aux__lbd*aux_flux*aux__fil, &
                                           lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
             lamb_effective = lamb_effective &
                            / IntegralALL( aux__lbd,aux__fil*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,IsOFFNum )
             deallocate( aux__fil )
          end if
          deallocate( aux__lbd )
          deallocate( aux_flux )
       end if
    end if

    return
END FUNCTION lamb_effective
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_lamb_effective( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_lamb_effective
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Thu Oct 13 12:53:07 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 012                                                          !
!  1) PropFilters
!  2) author_PropFilters
!  3) EvalTransmission
!  4) author_EvalTransmission
!  5) EvalTFluxes
!  6) author_EvalTFluxes
!  7) SLin_interp
!  8) author_SLin_interp
!  9) SLog_interp
! 10) author_SLog_interp
! 11) lamb_effective
! 12) author_lamb_effective

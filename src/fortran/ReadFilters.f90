! ###########################################################################
!     RESUME : Read Filters and retrieve ZERO point calibrations.           !
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
!                07) Numb_lbd  -> Number of points in filter                !
!                08) name_fil  -> Name of bandpass/filter                   !
!                09) Nfilters  -> Number of filters                         !
!                                                                           !
!     EXTRA ROUTINES : EvalTFluxes & IntegralALL                            !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadFilters( T_lambda,T_fluxes,T_l_Area,T_n_Area,magABsys,       &
                        magTGsys,standard,lamb_eff,Numb_lbd,name_fil,       &
                        file_dir,Nfilters,Int_Type,IsKeepOn,verbosity )
    use ModDataType
    implicit none
    real     (kind=RP), parameter :: cl_speed=299792.4580_RP
    real     (kind=RP), parameter :: fac4Pid2=1.19649518e40_RP
    
    integer  (kind=IB), intent(in) :: Int_Type
    integer  (kind=IB), intent(out) :: IsKeepOn,Nfilters

    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,i_filter,j_filter,NVecsize,N_lambda,     &
                          NVegalbd,NSun_lbd,NFBD_lbd,i_lambda,ilastnum,     &
                          ilastval
    integer  (kind=IB), dimension(:), intent(in out) :: Numb_lbd

    real     (kind=RP), target, dimension(:,:),intent(out) :: &
                                                                  T_lambda, &
                                                                  T_fluxes
    real     (kind=RP), pointer, dimension(:) :: P1lambda,P1fluxes
    real     (kind=RP), dimension(:), intent(in out) :: T_l_Area,           &
                                                        T_n_Area,           &
                                                        magABsys,           &
                                                        magTGsys,           &
                                                        standard,           &
                                                        lamb_eff

    real     (kind=RP), allocatable, dimension(:) :: lambveff,l_pivot2
    real     (kind=RP), allocatable, dimension(:) :: lambVega,fluxVega,     &
                                                     lamb_Sun,flux_Sun,     &
                                                     lamb_FBD,flux_FBD,     &
                                                     aux_flux,aux__lbd

    real     (kind=RP) :: lambda_i,lambda_f,aux_area,mAB_Vega,mTG_Vega,     &
                          flux_T01,flux_T02,flux_T03,lamb_T01,lamb_T02,     &
                          lamb_T03

    real     (kind=RP) :: delta_lbd
    
    character (len=*), intent(in) :: file_dir
    character (len=*), dimension(:), intent(in out) :: name_fil
    character (len=CH) :: axfilter,arq_lixo,arq_fil1,W1aux

    interface
        subroutine EvalTransmission( L_lambda,S_fluxes,T_lambda,T_fluxes,   &
                                     Ntlambda,fluxtran,N_lambda,IskeepOn,   &
                                     Int_Type,verbosity )
          use ModDataType
          implicit none
          integer  (kind=IB), intent(in out) :: IsKeepOn
          integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
          integer  (kind=IB), optional :: verbosity
          real     (kind=RP), dimension(:), intent(in) :: L_lambda,S_fluxes,&
                                                          T_lambda,T_fluxes
          real     (kind=RP), intent(in out) :: fluxtran
        end subroutine EvalTransmission
        subroutine EvalTFluxes( O_lambda,O_fluxes,T_lambda,T_fluxes,        &
                                T_l_area,N_lambda,fluxtran,lamb_eff,        &
                                Nspeclbd,IsKeepOn,verbosity )
            use ModDataType
            integer  (kind=IB), intent(in out) :: IsKeepOn
            integer  (kind=IB), intent(in) :: N_lambda,Nspeclbd
            integer  (kind=IB), optional :: verbosity
            real     (kind=RP), dimension(:), intent(in) :: O_lambda,       &
                                                            O_fluxes,       &
                                                            T_lambda,       &
                                                            T_fluxes
            real     (kind=RP), intent(in out) :: lamb_eff,fluxtran
            real     (kind=RP), intent(in) :: T_l_area
        end subroutine EvalTFluxes
        real (kind=RP) function IntegralALL( SXvalues,SYvalues,lambda_i,    &
                                             lambda_f,N_lambda,IsKeepOn,    &
                                             Int_Type,verbosity )
            use ModDataType
            integer  (kind=IB), intent(in out) :: IsKeepOn
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
           integer  (kind=IB), intent(in out) :: IsKeepOn
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

!  *** Reset Filter-STUFF ***************************************************
    T_l_Area = -999.0_RP
    T_n_Area = -999.0_RP
    magABsys = -999.0_RP
    magTGsys = -999.0_RP
    standard = -999.0_RP
    lamb_eff = -999.0_RP
    T_lambda = -999.0_RP
    T_fluxes = -999.0_RP

    name_fil = 'NONE'
!  *** Reset Filter-STUFF ***************************************************

    ilastnum = index(file_dir,' ') - 1

    NVecsize = size(T_lambda(:,1))

    allocate ( aux__lbd(NVecsize) )
    allocate ( aux_flux(NVecsize) )

! *** Read calibration stars ************************************************
!     RESUME : VEGA spectrum.                                               !
!              Intrinsic Flux: erg/s/cm2/A                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arq_fil1 = file_dir(1:ilastnum)//'VegaLR.dat'
    open  (21,status='old',file=arq_fil1,ERR=22)
    read  (21,*,ERR=22) arq_lixo,NVegalbd

    allocate ( lambVega(NVegalbd) )
    allocate ( fluxVega(NVegalbd) )

    do i_lambda=1,NVegalbd
        read  (21,*,ERR=22) lambVega(i_lambda),fluxVega(i_lambda)
    end do
    close (21)

! *** SUN spectrum **********************************************************
!     Intrinsic flux: erg/s/A.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arq_fil1 = file_dir(1:ilastnum)//'Sun_LR.dat'
    open(21,status='old',file=arq_fil1,ERR=22)
    read(21,*,ERR=22) arq_lixo,NSun_lbd

    allocate ( lamb_Sun(NSun_lbd) )
    allocate ( flux_Sun(NSun_lbd) )
    
    do i_lambda=1,NSun_lbd
        read(21,*,ERR=22) lamb_Sun(i_lambda),flux_Sun(i_lambda)
    end do
    close(21)

! *** Reading of the spectrum of BD+17o4708 *********************************
!     RESUME : F subdwarf used to calibrate the Thuan & Gunn system.        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !arq_fil1 = file_dir(1:ilastnum)//'BD+17o4708.dat'
    arq_fil1 = file_dir(1:ilastnum)//'BD+17d4708.dat'
    open(21,status='old',file=arq_fil1,ERR=22)
    read(21,*,ERR=22) arq_lixo,NFBD_lbd

    allocate ( lamb_FBD(NFBD_lbd) )
    allocate ( flux_FBD(NFBD_lbd) )
    
    do i_lambda=1,NFBD_lbd
        read(21,*,ERR=22) lamb_FBD(i_lambda),flux_FBD(i_lambda)
    end do
    close(21)
! *** Read calibration stars ************************************************

! *** Open ListFilters ******************************************************
    arq_fil1 = file_dir(1:ilastnum)//'ListFilters.txt'
    open  (unit=21,file=arq_fil1,status='old',ERR=23)
    read  (21,*,ERR=23) Nfilters

    if ( IsShowOn == 1 ) then
        write (*,'(4x,a)')   '[ReadFilters]'
        write (W1aux,'(i15)') Nfilters
        write (*,'(4x,a,a)') '... Nfilters: ',trim(adjustl(W1aux))
    end if

    do i_filter=1,Nfilters
        read  (21,*) axfilter
        name_fil(i_filter) = axfilter

        arq_fil1 = file_dir(1:ilastnum)//axfilter
        open  (unit=22,file=arq_fil1,status='old',ERR=23)
        do j_filter=1,15
            read  (22,*,END=20,ERR=23)
        end do
        do j_filter=1,NVecsize
           read  (22,*,END=20,ERR=23) T_lambda(j_filter,i_filter),          &
                                      T_fluxes(j_filter,i_filter)
        end do
20      Numb_lbd(i_filter) = j_filter-1
        close (22)
! *** Open ListFilters ******************************************************

    end do
    close (21)
! *** Open ListFilters ******************************************************

    !effective width
    !W = int(T dlamb) / max(T)
    
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
          
          lamb_eff(i_filter) = lambda_i
          standard(i_filter) = flux_T01
       else
          if (associated (P1lambda)) nullify( P1lambda )
          if (associated (P1fluxes)) nullify( P1fluxes )
          
          P1lambda => T_lambda(1:N_lambda,i_filter)
          P1fluxes => T_fluxes(1:N_lambda,i_filter)
          
          T_l_Area(i_filter) = IntegralALL( P1lambda,P1fluxes,lambda_i,     &
                                            lambda_f,N_lambda,IsKeepOn,     &
                                            Int_Type,00000000 )
          aux_area = T_l_Area(i_filter)
          
          aux_flux(1:N_lambda) = P1fluxes(1:N_lambda)                       &
                               / (P1lambda(1:N_lambda)                      &
                               * P1lambda(1:N_lambda)) * cl_speed
          aux_flux(1:N_lambda) = aux_flux(1:N_lambda) !* 1.0e13_RP
          
          T_n_Area(i_filter) = IntegralALL( P1lambda,aux_flux,lambda_i,     &
                                            lambda_f,N_lambda,IsKeepOn,     &
                                            Int_Type,00000000 )
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

          call EvalTransmission( lambVega,fluxVega,P1lambda,P1fluxes,       &
                                 N_lambda,flux_T01,NVegalbd,IskeepOn,       &
                                 Int_Type,0 )
          call EvalTransmission( lamb_Sun,flux_Sun,P1lambda,P1fluxes,       &
                                 N_lambda,flux_T02,NSun_lbd,IskeepOn,       &
                                 Int_Type,0 )
          call EvalTransmission( lamb_FBD,flux_FBD,P1lambda,P1fluxes,       &
                                 N_lambda,flux_T03,NFBD_lbd,IskeepOn,       &
                                 Int_Type,0 )

          !write (*,*) standard(i_filter),flux_T01
          standard(i_filter) = flux_T01
          lamb_eff(i_filter) = lamb_effective( P1lambda,P1fluxes,N_lambda,  &
                                               lambVega,fluxVega,NVegalbd,  &
                                               IskeepOn,Int_Type )

! *** Evaluate transmission flux ********************************************
       end if

! *** Magnitudes AB system **************************************************
       mAB_Vega = -2.5_RP                                                   &
                * log10( flux_T01*T_l_area(i_filter)/T_n_area(i_filter) )   &
                - 48.60_RP       
       magABsys(i_filter) = mAB_Vega       
! *** Magnitudes AB system **************************************************

! *** Magnitudes Thuan & Gunn system ****************************************
       if ( flux_T03 > 0.0_RP ) then
          mTG_Vega = -2.50_RP * log10(flux_T01/flux_T03) + 9.5_RP
       else
          mTG_Vega = -999._RP
       end if
       magTGsys(i_filter) = mTG_Vega
! *** Magnitudes Thuan & Gunn system ****************************************

! *** Print screen **********************************************************       
       if ( IsShowOn == 1_IB ) then          
          write (*,'(4x,a,a)') '... axfilter: ',trim(adjustl(name_fil(i_filter)))
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
          write (*,'(4x,a,a)') '++++++++++++++++++++++++++++++++++++++++++'
       end if
        
    end do
! *** Compute in Filters ****************************************************

    deallocate ( lambVega )
    deallocate ( fluxVega )
    deallocate ( lamb_Sun )
    deallocate ( flux_Sun )
    deallocate ( lamb_FBD )
    deallocate ( flux_FBD )
    deallocate ( aux__lbd )
    deallocate ( aux_flux )

    if ( IsShowOn == 1 ) then
        write (*,'(4x,a)')   '[ReadFilters]'
    end if

    return

! *** Issue warning - NO READING OF SPECTRUM ********************************
24  format(4x,65('?'))
22  write (*,24)
    write (*,'(4x,a)') '[ReadFilters] WARNING! NO SPECTRUM READ @@@@@@@'
    write (*,24)

23  write (*,24)
    write (*,'(4x,a)') '[ReadFilters] WARNING! PROBLEM - Filters List @'
    write (*,24)

    IsKeepOn = 0_IB
! *** Issue warning - NO EXTRAPOLATION **************************************

END SUBROUTINE ReadFilters
! ###########################################################################

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
!                09) verbosity -> Optional variable to print & check        !
!                                                                           !
!     OUTPUT   : 01) fluxtran -> Flux inside the filter                     !
!                02) lamb_eff -> Effective lambda                           !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvalTFluxes( O_lambda,O_fluxes,T_lambda,T_fluxes,T_l_area,       &
                        N_lambda,fluxtran,lamb_eff,Nspeclbd,IskeepOn,       &
                        verbosity )
    use ModDataType
    implicit none
    integer  (kind=IB), intent(in out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: N_lambda,Nspeclbd
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,jnum_aux,knum_aux
    real     (kind=RP), dimension(:), intent(in) :: O_lambda,O_fluxes,      &
                                                    T_lambda,T_fluxes
    real     (kind=RP), intent(in out) :: lamb_eff,fluxtran
    real     (kind=RP), intent(in) :: T_l_area
    real     (kind=RP) :: lambda_i,lambda_f,transinf,flux_inf,transsup,     &
                          flux_sup

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

! ###########################################################################
!     RESUME : Compute flux convoluted with the ransmission curve.          !
!                                                                           !
!     INPUT    : 01) O_lambda  -> Wavelength                                !
!                02) O_fluxes  -> Fluxes                                    !
!                03) T_lambda  -> Transmission lambda                       !
!                04) T_fluxes  -> Transmission                              !
!                05) T_l_area  -> Area computed with lambda                 !
!                06) N_lambda  -> Number of filter points                   !
!                07) Nspeclbd  -> Number of spectrum points                 !
!                08) IsKeepOn  -> Variable:       0 => Run aborted          !
!                09) verbosity -> Optional variable to print & check        !
!                                                                           !
!     OUTPUT   : 01) fluxtran  -> Flux inside the filter                    !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvalTransmission( L_lambda,S_fluxes,T_lambda,T_fluxes,Ntlambda,  &
                             fluxtran,N_lambda,IskeepOn,Int_Type,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,Nlbd_aux,ilastval,ind1,ind2

    real     (kind=RP), dimension(:), intent(in) :: L_lambda,S_fluxes,      &
                                                    T_lambda,T_fluxes
    real     (kind=RP), allocatable, dimension(:) :: aux__lbd,aux_flux,     &
                                                     aux__fil
    real     (kind=RP), intent(in out) :: fluxtran
    real     (kind=RP) :: lambda_i,lambda_f,dtlambda,d_lambda

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
         integer  (kind=IB), intent(in out) :: IsKeepOn
         integer  (kind=IB), intent(in) :: N_lambda,Int_Type
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,   &
                                                                SYvalues
         real     (kind=RP), intent(in) :: lambda_i,lambda_f
       end function IntegralALL
    end interface
    
    ! Interpolate to lowest dlambda
    dtlambda = minval( array=T_lambda(2:Ntlambda)-T_lambda(1:Ntlambda-1) )
    d_lambda = minval( array=L_lambda(2:N_lambda)-L_lambda(1:N_lambda-1),   &
                        mask=L_lambda(2:N_lambda) >= T_lambda(1) .AND.      &
                             L_lambda(2:N_lambda) <= T_lambda(Ntlambda) )

    !write (*,*) dtlambda,d_lambda
    
    ilastval = -999_IB
    if ( dtlambda < d_lambda ) then
       Nlbd_aux = Ntlambda
       allocate( aux__lbd(Ntlambda) )
       allocate( aux_flux(Ntlambda) )
       
       aux__lbd(1:Ntlambda) = T_lambda(1:Ntlambda)
       call LINinterpol( aux__lbd,aux_flux,Nlbd_aux,L_lambda,S_fluxes,N_lambda,ilastval,IsKeepOn,0 )
       
       lambda_i = T_lambda(1)
       lambda_f = T_lambda(Ntlambda)

       !flux = integrate[(f_lambda*lambda*T(lambda*(1+z))*dlambda)]
       !norm = integrate[(lambda*T(lambda*(1+z))*dlambda))]
       !photometry = flux/norm
        
       fluxtran = IntegralALL( aux__lbd,aux__lbd*T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
       fluxtran = fluxtran / IntegralALL( aux__lbd,aux__lbd*T_fluxes,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
       
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
             
             call LINinterpol( aux__lbd,aux__fil,Nlbd_aux,T_lambda,T_fluxes,Ntlambda,ilastval,IsKeepOn,0 )
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
             
             fluxtran = IntegralALL( aux__lbd,aux__lbd*aux_flux*aux__fil,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
             fluxtran = fluxtran / IntegralALL( aux__lbd,aux__lbd*aux__fil,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
             deallocate( aux__fil )
          end if
          deallocate( aux__lbd )
          deallocate( aux_flux )
       end if
    end if
    
    return
END SUBROUTINE EvalTransmission
! ###########################################################################

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
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2011                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SLin_interp( xorg_ini,xorg_fin,yorg_ini,yorg_fin,   &
                                     newvalue )
    use ModDataType
    implicit none
    real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,yorg_fin,  &
                                      newvalue

    SLin_interp = yorg_ini + (yorg_fin-yorg_ini) * (newvalue-xorg_ini)      &
                / (xorg_fin-xorg_ini)
    return
END FUNCTION SLin_interp
! ###########################################################################

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SLog_interp( xorg_ini,xorg_fin,yorg_ini,yorg_fin,   &
                                     newvalue )
    use ModDataType
    implicit none
    real     (kind=RP), parameter  :: setlimit=1.0e-37_RP
    real     (kind=RP), intent(in) :: xorg_ini,xorg_fin,yorg_ini,yorg_fin,  &
                                      newvalue

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

! ###########################################################################
!""" Unitwise Effective wavelength
!   leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
!   """
!   with Vega() as v:
!       s = self.reinterp(v.wavelength)
!       w = s._wavelength
!       leff = np.trapz(w * s.transmit * v.flux.magnitude, w, axis=-1)
!       leff /= np.trapz(s.transmit * v.flux.magnitude, w, axis=-1)
!   if self.wavelength_unit is not None:
!       return leff * unit[self.wavelength_unit]
!   else:
!       return leff
REAL (KIND=RP) FUNCTION lamb_effective( T_lambda,T_fluxes,Ntlambda,L_lambda,&
                                        S_fluxes,N_lambda,IskeepOn,Int_Type )
    use ModDataType
    implicit none
    integer  (kind=IB), intent(in out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: Ntlambda,N_lambda,Int_Type
    integer  (kind=IB) :: Nlbd_aux,ilastval,ind1,ind2
    
    real     (kind=RP), intent(in), dimension(Ntlambda) :: T_lambda,T_fluxes
    real     (kind=RP), intent(in), dimension(N_lambda) :: L_lambda,S_fluxes

    real     (kind=RP), allocatable, dimension(:) :: aux__lbd,aux_flux,aux__fil

    real     (kind=RP) :: dtlambda,d_lambda,lambda_i,lambda_f

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
       call LINinterpol( aux__lbd,aux_flux,Nlbd_aux,L_lambda,S_fluxes,N_lambda,ilastval,IsKeepOn,0 )
       
       lambda_i = T_lambda(1)
       lambda_f = T_lambda(Ntlambda)
       lamb_effective = IntegralALL( aux__lbd,aux__lbd*T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
       lamb_effective = lamb_effective / IntegralALL( aux__lbd,T_fluxes*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
       
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
             
             call LINinterpol( aux__lbd,aux__fil,Nlbd_aux,T_lambda,T_fluxes,Ntlambda,ilastval,IsKeepOn,0 )
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
             
             lamb_effective = IntegralALL( aux__lbd,aux__lbd*aux_flux*aux__fil,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
             lamb_effective = lamb_effective &
                            / IntegralALL( aux__lbd,aux__fil*aux_flux,lambda_i,lambda_f,Nlbd_aux,IsKeepOn,Int_Type,0 )
             deallocate( aux__fil )
          end if
          deallocate( aux__lbd )
          deallocate( aux_flux )
       end if
    end if

    return
END FUNCTION lamb_effective
! ###########################################################################

! Jean@Porto - Thu Oct 13 12:53:07 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!      PROGRAM GeneralTest
!        use ModDataType
!        implicit none
!        integer  (kind=IB), parameter :: Nb_max=1500,Nl_max=5000
!        integer  (kind=IB) :: IsKeepOn,Int_Type,iverbose,Nfilters
!        integer  (kind=IB), dimension(Nb_max) :: Numb_lbd
!        real     (kind=RP), dimension(Nl_max,Nb_max) :: T_lambda,T_fluxes
!        real     (kind=RP), dimension(Nb_max) :: T_l_Area,T_n_Area,magABsys,&
!                                                 magTGsys,standard,lamb_eff
!        character (len=CH), dimension(Nb_max) :: name_fil
!        character (len=CH) :: file_dir
!
!       interface
!          subroutine ReadFilters( T_lambda,T_fluxes,T_l_Area,T_n_Area,magABsys, &
!                                  magTGsys,standard,lamb_eff,Numb_lbd,name_fil, &
!                                  file_dir,Nfilters,Int_Type,IsKeepOn,verbosity )
!            use ModDataType
!            integer  (kind=IB), intent(in out) :: IsKeepOn,Nfilters
!            integer  (kind=IB), intent(in) :: Int_Type
!            integer  (kind=IB), optional :: verbosity
!            integer  (kind=IB), dimension(:), intent(in out) :: Numb_lbd
!            real     (kind=RP), target, dimension(:,:), intent(in out) ::       &
!                                                                      T_lambda, &
!                                                                      T_fluxes
!            real     (kind=RP), dimension(:), intent(in out) :: T_l_Area,       &
!                                                                T_n_Area,       &
!                                                                magABsys,       &
!                                                                magTGsys,       &
!                                                                standard,       &
!                                                                lamb_eff
!            character (len=*), intent(in) :: file_dir,file_out
!            character (len=*), dimension(:), intent(in out) :: name_fil
!          end subroutine ReadFilters
!       end interface
!
!       IsKeepOn = 1
!       iverbose = 1
!       file_dir = './'
!       file_out = 'Error.out'
!       Int_Type = 1
!       call ReadFilters( T_lambda,T_fluxes,T_l_Area,T_n_Area,magABsys, &
!                         magTGsys,standard,lamb_eff,Numb_lbd,name_fil, &
!                         file_dir,file_out,Nfilters,Int_Type,IsKeepOn, &
!                         iverbose )
!
!      END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 004                                                          !

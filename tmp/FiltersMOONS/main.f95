! A fortran95 program for G95
! By WQY
PROGRAM main

    !use ModDataType
    !use Modinterface
    !use ModVariableType
    use FAKEModules
    implicit none

    integer (kind=IB) :: inumeval,NVegalbd
    real    (kind=RP) :: L_Dis_pc
    real    (kind=RP), allocatable, dimension(:), target :: lambVega,fluxVega
    character(len=CH) :: arq_fil1

    logical, save :: first_entry = .true.

    IsKeepOn = 1                                  ! Flag to Keep On Going
    prolixon = 1                                  ! Print flag->Be prolix
    prolixof = 0                                  ! Not Print->Not prolix
    Read_Sun = 0                                  ! Read sun spectrum M2Lprolixon = 1

    gals_dir = '/run/media/jean/Isaac/FADO/D_FADO/REBETIKO/'

! *** Integration type ******************************************************
    Int_Type = 3
! *** Integration type ******************************************************

! *** Read filters and store on memory **************************************
    ilastnum = index(gals_dir,' ') - 1
    band_dir = gals_dir(1:ilastnum)//'F_FADO/'
    call ReadFilters( T_lambda,T_fluxes,T_l_Area,T_n_Area,magABsys,         &
                      magTGsys,standard,lamb_eff,Numb_lbd,name_fil,         &
                      band_dir,arqGlogs,Nfilters,Int_Type,IsKeepOn,         &
                      prolixon )

    !write (*,*) 'READ Filter'
! *** Read filters and store on memory **************************************

! *** Compute absolute magnitudes *******************************************
!     RESUME : Input spectrum must be in units of solar luminosities        !
!              per Ångström. The absolute magnitudes are stored in          !
!              MAGarray(Mslot_01,Mslot_02,Mslot_03), where:                 !
!                                                                           !
!              Mslot_01 -> Number of band-pass filters                      !
!              Mslot_02 -> Number of ZERO point calibrations                !
!              Mslot_03 -> Number of models                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! *** Read calibration stars ************************************************
!     RESUME : VEGA spectrum.                                               !
!              Intrinsic Flux: erg/s/cm2/A                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arq_fil1 = gals_dir(1:ilastnum)//'F_FADO/VegaLR.dat'
    open  (21,status='old',file=arq_fil1,ERR=22)
    read  (21,*,ERR=22) arq_lixo,NVegalbd

    allocate( lambVega(NVegalbd) )
    allocate( fluxVega(NVegalbd) )


    ! 7.68 ± 0.02 pc
    ! 1 pc = 3.08567758 × 10^(18)

    !L_Dis_pc = 7.68_RP !* 3.08567758e18
    L_Dis_pc = 10.0_RP
    log4PId2 = log10( 4.0_RP*Pi*L_Dis_pc**2 ) + 2.0_RP*log10( 3.086e18_RP )

    ! 40.12 ± 0.45[11] L☉
    do i_lambda=1,NVegalbd
        read  (21,*,ERR=22) lambVega(i_lambda),fluxVega(i_lambda)
        fluxVega(i_lambda) = fluxVega(i_lambda) * (10.0_RP**log4PId2) / (3.826e33_RP)
        !fluxVega(i_lambda) = fluxVega(i_lambda) !* (3.826e33_RP)
    end do
    close (21)

    write (*,*) log4PId2

! *** SUN spectrum **********************************************************
!     Intrinsic flux: erg/s/A.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arq_fil1 = gals_dir(1:ilastnum)//'F_FADO/Sun_LR.dat'
    open(21,status='old',file=arq_fil1,ERR=22)
    read(21,*,ERR=22) arq_lixo,NSun_lbd

    do i_lambda=1,NSun_lbd
        read(21,*,ERR=22) lamb_Sun(i_lambda),flux_Sun(i_lambda)
        ! --> Normalize to solar luminosity ---------------
        flux_Sun(i_lambda) = flux_Sun(i_lambda) / 3.826e33_RP
    end do
    close(21)

! Note that the AB system is defined such that a source with Fnu = 3.63 x 10-20 erg cm-2 s-1 Hz-1 has AB mag = 0 in every filter, and in general ABmag = - 2.5 log Fnu - 48.6.
!    flux_Sun(1:NSun_lbd) = 3.63e-20_RP * 2.99792458e18 / lamb_Sun(1:NSun_lbd)
    !L_Dis_pc = 10.0_RP
    !log4PId2 = log10( 4.0_RP*Pi*L_Dis_pc**2 ) + 2.0_RP*log10( 3.086e18_RP )
    !fluxVega(1:NSun_lbd) = 3.63e-20_RP * 2.99792458e18 / lambVega(1:NSun_lbd)**2 * (10.0_RP**log4PId2) / (3.826e33_RP) !* 1.7

    inumeval = 1

    if (first_entry) then
        nullify( MAG_spec ) ; nullify( P_lambda ) ; nullify( P1fluxes )
        first_entry = .false.
    end if

    if (associated (MAG_spec)) nullify( MAG_spec )
    MAG_spec => MAGarray(1:Np_max,1:Nf_max,inumeval)

    if (associated (P_lambda)) nullify( P_lambda )
    if (associated (P1fluxes)) nullify( P1fluxes )

    Int_Type = 3

    NGlambda = NVegalbd

    if ( .not. associated(P_lambda) ) then
        P_lambda => lambVega(1:NVegalbd)
    end if
    if ( .not. associated(P1fluxes) ) then
        P1fluxes => fluxVega(1:NVegalbd)
    end if
    !P_lambda => lambVega(1:NVegalbd)
    !P1fluxes => fluxVega(1:NVegalbd)

    call EvalFilters( P_lambda,P1fluxes,T_lambda,T_fluxes,T_l_Area,         &
                      T_n_Area,magABsys,magTGsys,standard,lamb_eff,         &
                      MAG_spec,Numb_lbd,name_fil,arq_logs,NGlambda,         &
                      Nfilters,Int_Type,IsKeepOn,prolixon )
! *** Compute absolute magnitudes *******************************************

    deallocate( lambVega )
    deallocate( fluxVega )

    !do i_lambda=1,7
    !    write (*,'(a10,7(e15.5))') name_fil(i_lambda),MAG_spec(i_lambda,:)
    !end do
    stop

24  format(4x,65('?'))
22  write (*,24)
    write (*,'(4x,a)') '[???????] WARNING! NO SPECTRUM READ @@@@@@@'
    write (*,24)

END PROGRAM main

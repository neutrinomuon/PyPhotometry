! ###########################################################################
!     RESUME : Integrate using different methods and spline1 interpolation  !
!              scheme. Original routine dates back to 2003-2004.            !
!                                                                           !
!              This is a Fortran 90 program for numerical integration of a  !
!              function using the Gaussian quadrature method. The input to  !
!              the function includes arrays of x-values (SXvalues) and      !
!              y-values (SYvalues) that define the function to be           !
!              integrated, a starting value (lambda_i) and an ending value  !
!              (lambda_f) for the integration interval, the number of       !
!              subintervals (N_lambda) to divide the interval into, an      !
!              integer (Int_Type) that defines the type of integration,     !
!              and an optional argument (verbosity) for the level of        !
!              detail to be output. The function returns a scalar value     !
!              (SIntegra) that represents the result of the integration.    !
!                                                                           !
!              The function makes use of two subroutines, LINinterpol and   !
!              GaussLegendreQuadrature. The LINinterpol subroutine          !
!              performs linear interpolation of a set of data, while the    !
!              GaussLegendreQuadrature subroutine performs the Gaussian     !
!              quadrature integration.                                      !
!                                                                           !
!              The Gaussian quadrature method is a numerical integration    !
!              technique that uses a set of weights and abscissas (points)  !
!              to approximate the definite integral of a function. The      !
!              method is considered to be very accurate, especially for     !
!              well-behaved functions.                                      !
!                                                                           !
!              The methods are given by Int_Type and may be summarized      !
!              bellow:                                                      !
!                                                                           !
!              ======================================================       !
!              |    Int_Type    |    Description                    |       !
!              |       0        |    R -> Right rectangle Integral  |       !
!              |       1        |    L -> Left rectangle Integral   |       !
!              |       2        |    T -> Trapezoidal rule          |       !
!              |       3        |    S -> Simple Integral           |       !
!              |       4        |    M -> Median rectangle Integral |       !
!              |       5        |    I -> Simpsonregel's rule       |       !
!              |       6        |    G -> Gauss-Legendre Quadrature |       !
!              ======================================================       !
!                                                                           !
!        OBS.: No extrapolation is done. Interpolation is only done         !
!              for the first and/or last points, depending if it is         !
!              included in the Sxvalues. The SPLINE1_Deg interpolation      !
!              was substituted by the LINEAR1_Deg subroutine, because       !
!              it is more efficient for a data set that is                  !
!              monotonically increasing.                                    !
!                                                                           !
!     INPUT  : 01) SXvalues -> 'X' (abcissas)  vector                       !
!              02) SYvalues -> 'Y' (ordenadas) vector                       !
!              03) lambda_i -> INITIAL 'X'                                  !
!              04) lambda_f -> FINAL   'X'                                  !
!              05) N_lambda -> # of elements in 'X' and 'Y'                 !
!              06) Int_Type -> Integration type                             !
!              07) Fverbose -> Optional variable to print & check           !
!                                                                           !
!     OUTPUT : 01) IntegralALL -> Area in a given interval                  !
!                                                                           !
!     LOG    : Error in Gaussian quadratute should be implemented.          !
!              Bug found in trapezium integration - May/14/2014             !
!                                                                           !
!     DEMO Tests : Very first tests                                         !
!                  Starting integration times - 10000 loops                 !
!                                                                           !
!     [IntegralALL]                                                         !
!     R-> Integral                                                          !
!            1       1     850    4985                                      !
!     ... SIntegra: 0.889702E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 2.1467 s                                                !
!     [EndingTimer]                                                         !
!     [IntegralALL]                                                         !
!     L-> Integral                                                          !
!            1       1     850    4985                                      !
!     ... SIntegra: 0.889445E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 2.1837 s                                                !
!     [EndingTimer]                                                         !
!     [IntegralALL]                                                         !
!     T-> Integral                                                          !
!            1       1     850    4985                                      !
!     ... SIntegra: 0.889573E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 2.6916 s                                                !
!     [EndingTimer]                                                         !
!     [IntegralALL]                                                         !
!     S-> Integral                                                          !
!            1       1     850    4985      4215.0000                       !
!     ... SIntegra: 0.893751E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 2.1487 s                                                !
!     [EndingTimer]                                                         !
!     [IntegralALL]                                                         !
!     M-> Integral                                                          !
!            1       1     850    4985                                      !
!     ... SIntegra: 0.889573E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 1.5988 s - Minimum integration time                     !
!     [EndingTimer]                                                         !
!     [IntegralALL]                                                         !
!     I-> Integral                                                          !
!            1       1     850    4985                                      !
!     ... SIntegra: 0.889573E+04                                            !
!     [IntegralALL]                                                         !
!     [EndingTimer]                                                         !
!     ... Telapsed: 1.7287 s                                                !
!     [EndingTimer]                                                         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Tue Nov 13 15:16:54  WET 2012                                !
!              Wed Jan  2 16:51:00  WET 2013                                !
!              Thu Jan  3 09:19:05  WET 2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION IntegralALL( SXvalues,SYvalues,lambda_i,lambda_f,   &
                                     N_lambda,IsKeepOn,Int_Type,verbosity )

    use ModDataType

    implicit none
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), intent(in) :: N_lambda,Int_Type
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: N_value1,N_value2,N_value3,N_value4,c_count1, &
                          c_count2,countdel,ilastval,IsShowOn
    integer  (kind=IB) :: IsOFF,IsONE,n_int
    real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,    &
                                                           SYvalues
    real     (kind=RP), dimension(1) :: aux1,aux2
    real     (kind=RP), allocatable, dimension(:) :: SLvalues,SRvalues, &
                                                     STvalues,SMvalues, &
                                                     SUvalues,SWvalues, &
                                                     SZvalues
    real     (kind=RP), intent(in) :: lambda_i,lambda_f
    real     (kind=RP) :: l_value1,fluxval1,l_value2,fluxval2,SIntegra, &
                          l_valmin,l_valmax,l_valaux,f_valaux,f_value1, &
                          f_value2,a,b
    character (len=CH) :: W1aux

    intrinsic adjustl, eoshift, count, maxval, min, minval, present,    &
              sum, trim
    
    interface
       subroutine LINinterpol( xx_value,yy_value,nxyvalue,xold_vec,     &
                               yold_vec,nold_vec,ilastval,IsKeepOn,     &
                               verbosity )
         use ModDataType
         integer  (kind=IB), intent(in out) :: ilastval
         integer  (kind=IB), intent(in) :: nold_vec, nxyvalue
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), dimension(0:nold_vec-1), intent(in) ::     &
                                                               xold_vec,&
                                                               yold_vec
         real     (kind=RP), dimension(0:nxyvalue-1), intent(in) ::     &
                                                               xx_value
         real     (kind=RP), dimension(0:nxyvalue-1), intent(out) ::    &
                                                               yy_value
       end subroutine LINinterpol
       subroutine GaussLegendreQuadrature( x,y,n_values,z,a,b,n_int )
         use ModDataType
         integer (kind=IB), intent(in) :: n_values
         integer (kind=IB), optional :: n_int
         real    (kind=RP), dimension(n_values), intent(in) :: x,y
         real    (kind=RP), intent(in)  :: a, b
         real    (kind=RP), intent(out) :: z
       end subroutine GaussLegendreQuadrature
    end interface

    !f2py real     (kind=RP), intent(in)  :: SXvalues,SYvalues
    !f2py integer  (kind=IB), intent(hide), depend(SXvalues) :: N_lambda=shape(SXvalues,0)

    !f2py real     (kind=RP), intent(in) :: lambda_i,lambda_f
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn

    ! Choose Trapezium as the default value
    !f2py integer  (kind=IB), optional :: Int_Type=2
    !f2py integer  (kind=IB), optional :: verbosity=0

    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    IsOFF = 0_IB
    IsONE = 1_IB
    
    SIntegra    = 0.0_RP                                                !Set integral to 0 @@@@@@@@@@@@@@@@
    IntegralALL = SIntegra                                              !Set integral to 0 @@@@@@@@@@@@@@@@

    l_valmin = minval( SXvalues(1:N_lambda) )
    l_valmax = maxval( SXvalues(1:N_lambda) )

    allocate( SLvalues(N_lambda) )
    allocate( SRvalues(N_lambda) )
    allocate( STvalues(N_lambda) )
    allocate( SMvalues(N_lambda) )
    allocate( SUvalues(N_lambda) )
    allocate( SWvalues(N_lambda) )
    allocate( SZvalues(N_lambda) )

    SLvalues = 0.0_RP
    SRvalues = 0.0_RP
    STvalues = 0.0_RP
    SMvalues = 0.0_RP
    SUvalues = 0.0_RP
    SWvalues = 0.0_RP
    SZvalues = 0.0_RP

    if ( lambda_i >= l_valmax .OR. lambda_f <= l_valmin ) then
        if ( IsShowOn == 1_IB ) then 
           write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
           write (*,'(4x,a)')  '[IntegralALL] l & F is out of range @@'
        end if
        SIntegra    = -999.0_RP
        IntegralALL = SIntegra
        IsKeepOn    = 0_IB
        return
    end if

    countdel = count( SXvalues(1:N_lambda-1) >= SXvalues(2:N_lambda) )
    if ( countdel > 0 ) then
        if ( IsShowOn == 1_IB ) then 
           write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
           write (*,'(4x,a)')  '[IntegralALL] l is not in increasing @'
        end if
        SIntegra    = -999.0_RP
        IntegralALL = SIntegra
        IsKeepOn    = 0_IB
        return
    end if

! *** Issue warnings - NO EXTRAPOLATION *************************************
10  format(4x,65('?'))
    if ( lambda_i < l_valmin ) then
        if ( IsShowOn == 1_IB ) then 
           write (*,10)
           write (*,'(4x,a)') '[IntegralALL] WARNING! lambda_i is out of range'
           write (*,'(4x,a)') '[IntegralALL] WARNING! NO EXTRAPOLATION @@@@@@@'
           write (*,'(4x,2(a,f15.3))') '[IntegralALL] ',lambda_i,' < ',l_valmin
           write (*,10)
        end if
    end if
    if ( lambda_f > l_valmax ) then
        if ( IsShowOn == 1_IB ) then 
           write (*,10)
           write (*,'(4x,a)') '[IntegralALL] WARNING! lambda_f is out of range'
           write (*,'(4x,a)') '[IntegralALL] WARNING! NO EXTRAPOLATION @@@@@@@'
           write (*,'(4x,2(a,f15.3))') '[IntegralALL] ',lambda_f,' > ',l_valmax
           write (*,10)
        end if
    end if
! *** Issue warnings - NO EXTRAPOLATION *************************************

    if ( IsShowOn == 1_IB ) then
        write (*,'(4x,a)') '[IntegralALL]'
    end if

    select case (Int_Type)
    case default
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'S-> Integral'
        end if

        c_count1 = 1_IB
        c_count2 = 1_IB
        do while ( c_count1 <= N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1) <= lambda_f ) then
                c_count2 = c_count1
                if ( c_count1 < N_lambda ) then
                    STvalues(c_count1) = SXvalues(c_count1+1)-SXvalues(c_count1)
                else
                    STvalues(c_count2) = STvalues(c_count2-1)
                end if
            end if
            c_count1 = c_count1 + 1
        end do

        !SIntegra = sum(SYvalues*STvalues, SXvalues.GE.lambda_i.AND.SXvalues.LE.lambda_f)

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )

        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8),f15.4)') N_value1,N_value2,N_value3,N_value4,SXvalues(c_count2)
        end if

        if ( N_value2 == 1_IB ) then
            STvalues(c_count2) = STvalues(c_count2-1)
        end if

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i

            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value1 = aux2(1)
            l_valaux = minval( SXvalues, SXvalues > l_value1 )

            SIntegra = SIntegra + f_value1 * (l_valaux - l_value1)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,f_value1,(f_value1*(l_valaux-l_value1))
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            STvalues(c_count2) = 0.0_RP
            l_value2 = lambda_f

            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)
            l_valaux = maxval( SXvalues, SXvalues < l_value2 )

            STvalues(c_count2) = (l_value2 - l_valaux)

            SIntegra = SIntegra + f_value2 * (l_value2 - l_valaux)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,f_value2,(f_value2*(l_value2-l_valaux))
            end if
        end if

        SIntegra = SIntegra + sum(SYvalues * STvalues)

    case (0)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'R-> Integral'
        end if

        c_count1 = 1_IB
        c_count2 = 1_IB
        do while ( c_count1 < N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1) <= lambda_f ) then
                SRvalues(c_count1) = SXvalues(c_count1+1)-SXvalues(c_count1)
                c_count2 = c_count1
            end if
            c_count1 = c_count1 + 1
        end do

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8))') N_value1,N_value2,N_value3,N_value4
        end if

        if ( N_value2 == 1_IB .AND. SXvalues(N_lambda) == lambda_f ) then
            SRvalues(c_count2+1) = 0.0_RP
        else
            if ( SXvalues(N_lambda) < lambda_f ) then
                SRvalues(c_count2+1) = 0.0_RP
            else
                SRvalues(c_count2)   = 0.0_RP
            end if
        end if

        !write (*,*) c_count2,c_count1,SXvalues(c_count2+1),SXvalues(c_count2), &
        !            SRvalues(c_count2),SRvalues(c_count2+1),lambda_i,lambda_f

        SIntegra = sum( SYvalues*SRvalues )

        !write (*,*) SIntegra

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i

            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            fluxval1 = aux2(1)
            
            l_valaux = minval( SXvalues, SXvalues > l_value1 )

            SIntegra = SIntegra + (l_valaux-l_value1) * fluxval1
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,fluxval1, &
                                                      ((l_valaux-l_value1)*fluxval1)
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            l_value2 = lambda_f

            l_valaux = maxval( SXvalues, SXvalues < l_value2 )
            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_valaux = aux2(1)

            SIntegra = SIntegra + (l_value2-l_valaux) * f_valaux
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value1,fluxval1, &
                                                      ((l_value2-l_valaux)*f_valaux)
            end if
        end if

    case (1)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'L-> Integral'
        end if

        c_count1 = 2_IB
        c_count2 = 1_IB
        do while ( c_count1 <= N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1) <= lambda_f ) then
                SLvalues(c_count1) = SXvalues(c_count1) - SXvalues(c_count1-1)
            end if
            if ( SXvalues(c_count1) <= lambda_i ) then
                c_count2 = c_count1
            end if
            c_count1 = c_count1 + 1
        end do

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8))') N_value1,N_value2,N_value3,N_value4
        end if

        if ( N_value1 == 1_IB ) then
            SLvalues(c_count2+0) = 0.0_RP
        end if
        if ( N_value1 /= 1_IB ) then
            SLvalues(c_count2+1) = 0.0_RP
        end if
        !write (*,*) c_count2,c_count1,SXvalues(c_count2),SXvalues(c_count2-1),SLvalues(c_count2),SLvalues(c_count2-1),lambda_i,lambda_f

        SIntegra = sum( SYvalues*SLvalues )

        !write (*,*) SIntegra

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i

            l_valaux = minval( SXvalues, SXvalues > l_value1 )
            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_valaux = aux2(1)

            SIntegra = SIntegra + (l_valaux-l_value1) * f_valaux
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,fluxval1,((l_valaux-l_value1)*f_valaux)
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            l_value2 = lambda_f
            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)

            l_valaux = maxval( SXvalues, SXvalues < l_value2 )

            SIntegra = SIntegra + (l_value2-l_valaux) * f_value2
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,f_value2,((l_value2-l_valaux)*f_value2)
            end if
        end if

    case (2)
        if ( IsShowOn == 1_IB) then
            write (*,'(4x,a)') 'T-> Integral'
        end if

        c_count1 = 1_IB
        c_count2 = 1_IB
        do while ( c_count1 < N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1) <= lambda_f ) then
                STvalues(c_count1) = SXvalues(c_count1+1)-SXvalues(c_count1)
                c_count2 = c_count1
            end if
            c_count1 = c_count1 + 1
        end do

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8))') N_value1,N_value2,N_value3,N_value4
        end if

        if ( N_value2 == 1_IB .AND. SXvalues(N_lambda) == lambda_f ) then
            STvalues(c_count2+1) = 0.0_RP
        else
            if ( SXvalues(N_lambda) < lambda_f ) then
                STvalues(c_count2+1) = 0.0_RP
            else
                STvalues(c_count2)   = 0.0_RP
            end if
        end if

        SWvalues = eoshift( SYvalues, shift=1, boundary=0.0_RP )
        SZvalues(:) = SWvalues(:) + SYvalues(:)

        !write (*,*) c_count2,c_count1,SXvalues(c_count2+1),SXvalues(c_count2),SRvalues(c_count2+1),lambda_i,lambda_f

        SIntegra = sum( 0.5_RP*SZvalues*STvalues )

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i
            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value1 = aux2(1)

            l_valaux = minval( SXvalues, SXvalues > l_value1 )
            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_valaux = aux2(1)

            SIntegra = SIntegra + 0.5_RP * (f_value1 + f_valaux) * (l_valaux - l_value1)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,&
                             fluxval1,(0.5_RP*(f_value1 + f_valaux)*(l_valaux-l_value1))
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            l_value2 = lambda_f
            l_valaux = maxval( SXvalues, SXvalues < l_value2 )
            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_valaux = aux2(1)

            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)

            SIntegra = SIntegra + 0.5_RP * (f_value2 + f_valaux) * (l_value2 - l_valaux)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,&
                               f_value2,(0.5_RP*(f_value2+f_valaux)*(l_value2-l_valaux))
            end if
        end if

    case (3)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'S-> Integral'
        end if

        c_count1 = 1
        c_count2 = 1
        do while ( c_count1 <= N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1) <= lambda_f ) then
                c_count2 = c_count1
                if ( c_count1 < N_lambda ) then
                    STvalues(c_count1) = SXvalues(c_count1+1)-SXvalues(c_count1)
                 else
                    STvalues(c_count2) = STvalues(c_count2-1)
                 end if
            end if
            c_count1 = c_count1 + 1
        end do

        !SIntegra = sum(SYvalues*STvalues, SXvalues.GE.lambda_i.AND.SXvalues.LE.lambda_f)

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8),f15.4)') N_value1,N_value2,N_value3,N_value4,SXvalues(c_count2)
        end if

        if ( N_value2 == 1_IB ) then
            STvalues(c_count2) = STvalues(c_count2-1)
        end if

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i

            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value1 = aux2(1)

            l_valaux = minval( SXvalues, SXvalues > l_value1 )

            SIntegra = SIntegra + f_value1 * (l_valaux - l_value1)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,f_value1,(f_value1*(l_valaux-l_value1))
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            STvalues(c_count2) = 0.0_RP
            l_value2 = lambda_f

            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)

            l_valaux = maxval( SXvalues, SXvalues < l_value2 )

            STvalues(c_count2) = (l_value2 - l_valaux)

            SIntegra = SIntegra + f_value2 * (l_value2 - l_valaux)
            if ( IsShowOn == 1 ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,f_value2,(f_value2*(l_value2-l_valaux))
            end if
        end if

        SIntegra = SIntegra + sum( SYvalues*STvalues )

    case (4)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'M-> Integral'
        end if

        c_count1 = 1
        c_count2 = 1
        ilastval = -999 ! *** Setting interpolation to -999 @@@@
        do while ( c_count1 < N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1+1) <= lambda_f ) then
                SMvalues(c_count1) = SXvalues(c_count1+1)-SXvalues(c_count1)
                SWvalues(c_count1) = SXvalues(c_count1) + SMvalues(c_count1) * 0.5_RP
                l_valaux = SWvalues(c_count1)

                aux1(1) = l_valaux
                call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
                f_valaux = aux2(1)

                SZvalues(c_count1) = f_valaux

                c_count2 = c_count1
            end if
            c_count1 = c_count1 + 1
        end do

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1 ) then
            write (*,'(4x,4(i8))') N_value1,N_value2,N_value3,N_value4
        end if
        !write (*,*) c_count2,c_count1,SXvalues(c_count2),SXvalues(c_count2-1),SMvalues(c_count2),SMvalues(c_count2-1),lambda_i,lambda_f

        SIntegra = sum( SMvalues*SZvalues )

        !write (*,*) SIntegra

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i
            l_valaux = minval( SXvalues, SXvalues > l_value1 )

            l_value2 = l_value1 + (l_valaux - l_value1) * 0.5_RP
            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)

            SIntegra = SIntegra + f_value2 * (l_valaux-l_value1)
            if ( IsShowOn == 1_IB ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,f_value2, &
                                                       f_value2*(l_valaux-l_value1)
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            l_value2 = lambda_f
            l_valaux = maxval( SXvalues, SXvalues < l_value2 )

            l_value1 = l_valaux + (l_value2-l_valaux) * 0.5_RP
            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value1 = aux2(1)

            SIntegra = SIntegra + f_value1 * (l_value2-l_valaux)
            if ( IsShowOn == 1 ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,f_value1, &
                                                       f_value1*(l_value2-l_valaux)
            end if
        end if

    case (5)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'I-> Integral'
        end if

        c_count1 = 1_IB
        c_count2 = 1_IB
        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        do while ( c_count1 < N_lambda )
            if ( SXvalues(c_count1) >= lambda_i .AND. SXvalues(c_count1+1) <= lambda_f ) then
                SMvalues(c_count1) = (SXvalues(c_count1+1)-SXvalues(c_count1))/6.0_RP
                SWvalues(c_count1) = SYvalues(c_count1+1)
                l_valaux = (SXvalues(c_count1+1)+SXvalues(c_count1))*0.5_RP

                aux1(1) = l_valaux
                call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
                f_valaux = aux2(1)

                SZvalues(c_count1) = f_valaux

                c_count2 = c_count1
            end if
            c_count1 = c_count1 + 1
        end do

        N_value1 = count( SXvalues == lambda_i )
        N_value2 = count( SXvalues == lambda_f )
        N_value3 = count( SXvalues  < lambda_i )
        N_value4 = count( SXvalues  > lambda_f )
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,4(i8))') N_value1,N_value2,N_value3,N_value4
        end if
        !write (*,*) c_count2,c_count1,SXvalues(c_count2),SXvalues(c_count2-1),SMvalues(c_count2),SMvalues(c_count2-1),lambda_i,lambda_f

        SIntegra = sum( SMvalues*(SYvalues+SWvalues+4.0_RP*SZvalues) )

        !write (*,*) SIntegra

        ilastval = -999_IB ! *** Setting interpolation to -999 @@@@
        if ( N_value1 == 0_IB .AND. N_value3 > 0_IB ) then
            l_value1 = lambda_i
            l_valaux = minval( SXvalues, SXvalues > l_value1 )

            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            fluxval1 = aux2(1)

            l_value2 = (l_value1+l_valaux)*0.5_RP
            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value2 = aux2(1)

            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            fluxval2 = aux2(1)

            SIntegra = SIntegra + (l_valaux-l_value1)/6.0_RP*(fluxval1+fluxval2+4.0_RP*f_value2)
            if ( IsShowOn == 1 ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 1',l_valaux,l_value1,f_value2, &
                                  (l_valaux-l_value1)/6.0_RP*(fluxval1+fluxval2+4.0_RP*f_value2)
            end if
        end if

        if ( N_value2 == 0_IB .AND. N_value4 > 0_IB ) then
            l_value2 = lambda_f
            l_valaux = maxval( SXvalues, SXvalues < l_value2 )

            aux1(1) = l_valaux
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            fluxval2 = aux2(1)

            l_value1 = (l_value2+l_valaux)*0.5_RP
            aux1(1) = l_value1
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            f_value1 = aux2(1)

            aux1(1) = l_value2
            call LINinterpol( aux1,aux2,IsONE,SXvalues,SYvalues,N_lambda,ilastval,IsKeepOn,IsOFF )
            fluxval1 = aux2(1)

            SIntegra = SIntegra + (l_value2-l_valaux)/6.0_RP*(fluxval1+fluxval2+4.0_RP*f_value1)
            if ( IsShowOn == 1 ) then
                write (*,'(4x,a,4(f15.4))') '++++++ 2',l_valaux,l_value2,f_value1, &
                                  (l_valaux-l_value2)/6.0_RP*(fluxval1+fluxval2+4.0_RP*f_value1)
            end if
        end if

    case (6)
        if ( IsShowOn == 1_IB ) then
            write (*,'(4x,a)') 'G-> Integral'
        end if

        a = max(lambda_i,SXvalues(1))
        b = min(lambda_f,SXvalues(N_lambda))
        n_int = 20
        call GaussLegendreQuadrature( SXvalues,SYvalues,N_lambda,SIntegra,a,b,n_int )
    end select

    if ( IsShowOn == 1_IB ) then
        write (W1aux,'(e25.6)') SIntegra
        write (*,'(4x,a,a)') '... SIntegra: ',trim(adjustl(W1aux))
        write (*,'(4x,a)') '[IntegralALL]'
    end if

    deallocate ( SLvalues )
    deallocate ( SRvalues )
    deallocate ( STvalues )
    deallocate ( SMvalues )
    deallocate ( SUvalues )
    deallocate ( SWvalues )
    deallocate ( SZvalues )

    IntegralALL = SIntegra
    return
END FUNCTION IntegralALL
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_IntegralALL( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_IntegralALL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Wed Sep 28 11:30:26 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!     PROGRAM GeneralTest
!       use ModDataType
!       implicit none
!       integer  (kind=IB), parameter :: Nl_max=80000
!       integer  (kind=IB) :: IsKeepOn,Int_Type,Intloops,i
!       real     (kind=RP), dimension(Nl_max) :: l,f,e,m
!       real     (kind=RP) :: int,lambda_i,lambda_f,ll1,ff1,ll2,ff2,ll3,ff3,start
!       character (len=CH) :: arq_gals
!
!       interface
!          subroutine CleanAClear( a,j,k )
!            use ModDataType
!            integer  (kind=IB), optional :: j,k
!            real     (kind=RP), dimension(:), intent (in out) :: a
!          end subroutine CleanAClear
!          real (kind=RP) function IntegralALL( SXvalues,SYvalues, &
!               lambda_i,lambda_f, &
!               N_lambda,IsKeepOn, &
!               Int_Type,verbosity )
!            use ModDataType
!            integer  (kind=IB), intent(in out) :: IsKeepOn
!            integer  (kind=IB), intent(in) :: N_lambda,Int_Type
!            integer  (kind=IB), optional :: verbosity
!            real     (kind=RP), dimension(N_lambda), intent(in) :: SXvalues,SYvalues
!            real     (kind=RP), intent(in) :: lambda_i,lambda_f
!          end function IntegralALL
!          subroutine Start_Timer( start,k )
!            use ModDataType
!            integer  (kind=IB), optional :: k
!            real     (kind=RP), intent(in out) :: start
!          end subroutine Start_Timer
!          subroutine EndingTimer( start,k )
!            use ModDataType
!            integer  (kind=IB), optional :: k
!            real     (kind=RP), intent(in out) :: start
!          end subroutine EndingTimer
!       end interface
!
!       l = 0.0_RP
!       f = 0.0_RP
!       e = 0.0_RP
!       m = 0.0_RP
!
!       !lambda_i = 3645.000_RP
!       !lambda_f = 8859.000_RP
!
!       l(1) = 4000.000_RP
!       l(2) = 4100.000_RP
!       l(3) = 4150.000_RP
!       l(4) = 4200.000_RP
!       l(5) = 4230.000_RP
!
!       f(1) = 010.000_RP
!       f(2) = 005.000_RP
!       f(3) = 008.000_RP
!       f(4) = 006.100_RP
!       f(5) = 007.500_RP
!
!       ll1   = 4050.000_RP
!       ff1   = f(1) + (f(2)-f(1)) / (l(2)-l(1)) * (ll1 - l(1))
!
!       ll2   = 4215.000_RP
!       ff2   = f(4) + (f(5)-f(4)) / (l(5)-l(4)) * (ll2 - l(4))
!
!       !R
!       int = (l(2)-l(1))*f(1) + (l(3)-l(2))*f(2) + (l(4)-l(3))*f(3) !+ (l(5)-l(4))*f(4)
!       int =  ff1 * (l(2) - ll1) +  (l(3)-l(2))*f(2) + (l(4)-l(3))*f(3) + (ll2-l(4))*f(4)               !+ (l(5)-l(4))*f(4)
!       write (*,'(4x,a,f15.5)') 'R-> ',int
!
!       !L
!       int = (l(2)-l(1))*f(2) + (l(3)-l(2))*f(3) + (l(4)-l(3))*f(4) !+ (l(5)-l(4))*f(5)
!       int = (l(2)-ll1)*f(2) + (l(3)-l(2))*f(3) + (l(4)-l(3))*f(4) + (ll2-l(4))*ff2       !+ (l(5)-l(4))*f(5)
!       write (*,'(4x,a,f15.5)') 'L-> ', int
!
!       !T
!       int = 0.50_RP * ((l(2)-l(1))*(f(1)+f(2)) + (l(3)-l(2))*(f(2)+f(3)) + (l(4)-l(3))*(f(3)+f(4))) !&
!       !+ 0.5_RP * ((l(5)-l(4))*(f(4)+f(5)))
!       int = 0.50_RP * ((l(3)-l(2))*(f(2)+f(3)) + (l(4)-l(3))*(f(3)+f(4))) + 0.50_RP * ((l(2)-ll1)*(f(2)+ff1)) &
!            +0.50_RP * (ll2-l(4))*(f(4)+ff2)
!       write (*,'(4x,a,f15.5)') 'T-> ', int
!
!       !S
!       int = (l(2)-l(1))*f(1) + (l(3)-l(2))*f(2) + (l(4)-l(3))*f(3) + (l(4)-l(3))*f(4) !+ (l(5)-l(4))*f(5) !+ (l(5)-l(4))*f(5)
!       int = (l(2)-ll1)*ff1 + (l(3)-l(2))*f(2) + (l(4)-l(3))*f(3) + (ll2-l(4))*f(4) + (ll2-l(4))*ff2          !+ (l(5)-l(4))*f(5) !+ (l(5)-l(4))*f(5)
!       write (*,'(4x,a,f15.5)') 'S-> ',int
!
!       !M
!       !ll3 = l(1) + (l(2)-l(1)) / 2.0_RP
!       !ff3 = f(1) + (f(2)-f(1)) / (l(2)-l(1)) * (ll3 - l(1))
!       !int = ff3 * (l(2)-l(1))
!
!       ll3 = ll1 + (l(2)-ll1) / 2.0_RP
!       ff3 = ff1 + (f(2)-ff1) / (l(2)-ll1) * (ll3 - ll1)
!       int = ff3 * (l(2)-ll1)
!
!       ll3 = l(2) + (l(3)-l(2)) / 2.0_RP
!       ff3 = f(2) + (f(3)-f(2)) / (l(3)-l(2)) * (ll3 - l(2))
!       int = int + ff3 * (l(3)-l(2))
!
!       ll3 = l(3) + (l(4)-l(3)) / 2.0_RP
!       ff3 = f(3) + (f(4)-f(3)) / (l(4)-l(3)) * (ll3 - l(3))
!       int = int + ff3 * (l(4)-l(3))
!
!       !ll3 = l(4) + (l(5)-l(4)) / 2.0_RP
!       !ff3 = f(4) + (f(5)-f(4)) / (l(5)-l(4)) * (ll3 - l(4))
!       !int = int + ff3 * (l(5)-l(4))
!
!       ll3 = l(4) + (ll2-l(4)) / 2.0_RP
!       ff3 = f(4) + (ff2-f(4)) / (ll2-l(4)) * (ll3 - l(4))
!       int = int + ff3 * (ll2-l(4))
!       write (*,'(4x,a,f15.5)') 'M-> ',int
!
!       !I
!       ll3 = (l(2)+ll1) / 2.0_RP
!       ff3 = ff1 + (f(2)-ff1) / (l(2)-ll1) * (ll3 - ll1)
!       int = (l(2)-ll1)/6.0_RP*(ff1+f(2)+4.0_RP*ff3)
!
!       ll3 = (l(3)+l(2)) / 2.0_RP
!       ff3 = f(2) + (f(3)-f(2)) / (l(3)-l(2)) * (ll3 - l(2))
!       int = int + (l(3)-l(2))/6.0_RP*(f(2)+f(3)+4.0_RP*ff3)
!
!       ll3 = (l(4)+l(3)) / 2.0_RP
!       ff3 = f(3) + (f(4)-f(3)) / (l(4)-l(3)) * (ll3 - l(3))
!       int = int + (l(4)-l(3))/6.0_RP*(f(3)+f(4)+4.0_RP*ff3)
!
!       ll3 = (ll2+l(4)) / 2.0_RP
!       ff3 = f(4) + (ff2-f(4)) / (ll2-l(4)) * (ll3 - l(4))
!       int = int + (ll2-l(4))/6.0_RP*(f(4)+ff2+4.0_RP*ff3)
!
!       write (*,'(4x,a,f15.5)') 'I-> ',int
!
!       i = 5
!       lambda_i = 4050.000_RP
!       lambda_f = 4215.000_RP
!
!       IsKeepOn = 1
!       Int_Type = 0
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 1
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 2
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 3
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 4
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 5
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       Int_Type = 6
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!
!       write (*,*)
!       write (*,'(4x,a)') 'Starting integration times - 10000 loops'
!
!       l = 0.0_RP
!       f = 0.0_RP
!       e = 0.0_RP
!       m = 0.0_RP
!
!       arq_gals = '/home/jean/FADO/C_FADO/Input/0266.51602.089.7xt'
!       open  (unit=1,file=arq_gals,status='old')
!       do i=1,size(l)
!          read (1,*,END=20,ERR=20) l(i),f(i)
!       end do
!       i = i - 1
!       close (1)
!
!       do Int_Type=0,6
!          call Start_Timer( start,0 )
!          do Intloops=1,10000
!             int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,0)
!          end do
!          int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!          call EndingTimer( start,1 )
!       end do
!
!       write (*,*)
!       l(2) = l(3)
!       int = IntegralALL(l,f,lambda_i,lambda_f,i,IsKeepOn,Int_Type,1)
!       write (*,'(4x,a,f15.5)') 'P-> ',int
!       !write (*,*) (lambda_f**4/4.0_RP - lambda_i**4/4.0_RP),int-(lambda_f**4/4.0_RP - lambda_i**4/4.0_RP)
!
!     END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 002                                                          !
!
! 1) IntegralALL
! 2) author_IntegralALL


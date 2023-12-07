! ###########################################################################
!     RESUME : Numerical Gauss-Legendre Quadrature. The subroutine          !
!              calculates an approximation of the definite integral of      !
!              the function represented by x and y over the interval        !
!              [a,b] using the Gauss-Legendre quadrature method. The        !
!              result is returned in the z output variable.                 !
!                                                                           !
!              The subroutine uses the LINinterpol subroutine to            !
!              perform linear interpolation of the function values at       !
!              new points required for the Gaussian quadrature              !
!              calculation. The subroutine calculates the Gaussian          !
!              quadrature points and weights using a recursive formula      !
!              for the Legendre polynomials. Finally, the integral          !
!              approximation is calculated as the sum of the product        !
!              of the weights and the function values at the                !
!              quadrature points, multiplied by (b-a)/2.                    !
!                                                                           !
!     INPUT  : 01) x        -> 'X' (abcissas)  vector                       !
!              02) y        -> 'Y' (ordenadas) vector                       !
!              03) n_values -> # of elements in 'X' and 'Y'                 !
!              04) a        -> INITIAL 'X'                                  !
!              05) b        -> FINAL   'X'                                  !
!              06) n_int    -> Optional variable with number of iterations  !
!                                                                           !
!     OUTPUT : 01) z -> Area in a given interval using Gauss-Legendre       !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Tue Nov 13 15:16:54  WET 2012                                !
!              Wed Jan  2 16:51:00  WET 2013                                !
!              Thu Jan  3 09:19:05  WET 2013                                !
!                                                                           !
!         LOG: Ter  5 Dez 2023 11:59:26 WET                                 !
!              Problem with definition of variable n_interactions           !
!              Corrected the problem by changing n_interactions to n_int    !
!              if ( n_int > 1_IB ) then                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GaussLegendreQuadrature( x,y,n_values,z,a,b,n_int )

    use ModDataType

    implicit none
    integer (kind=IB)             :: n = 10, k, npt, ilastval,IsKeepOn,     &
                                     IsOFF, n_interations
    integer (kind=IB), intent(in) :: n_values
    integer (kind=IB), optional :: n_int

    real    (kind=RP), dimension(n_values), intent (in) :: x,y
    real    (kind=RP), intent (in)  :: a, b
    real    (kind=RP), intent (out) :: z
    
    real    (kind=RP), allocatable, dimension(:,:) :: r
    real    (kind=RP), allocatable, dimension(:) :: y_alloc,x_alloc

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
    end interface

    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py integer  (kind=IB), intent(hide), depend(x) :: n_values=shape(x,0)

    !f2py real     (kind=RP), intent(in) :: a,b
    !f2py real     (kind=RP), intent(out) :: z

    !f2py integer  (kind=IB), optional :: n_int=20
    
    IsKeepOn = 0_IB
    IsOFF = 0_IB

    n_interations = 20_IB
    if ( present(n_int) ) then
       if ( n_int > 1_IB ) then
          n_interations = n_int
       else
          n_interations = 20_IB
       end if
    else
       n_interations = 20_IB
    end if

    do n=1,n_int
       r = gaussquadrature(n)
       
       npt = size(r(1,:))
       allocate( x_alloc(npt) )
       allocate( y_alloc(npt) )
       
       ! Here interpolate to the new values required
       ilastval = -999 ! *** Setting interpolation to -999 @@@@
       x_alloc = (a+b)/2+r(1,:)*(b-a)/2
       call LINinterpol( x_alloc,y_alloc,npt,x,y,n_values,ilastval,IsKeepOn,IsOFF )
       
       z = (b-a)/2 * dot_product( r(2,:),y_alloc )

       deallocate( r )
       deallocate( x_alloc )
       deallocate( y_alloc )
    end do
         
  contains 

    FUNCTION gaussquadrature(n) result(r)
        use ModDataType  
        real    (kind=RP), parameter :: pi = 4.0_RP*atan(1.0_RP)
        
        integer (kind=IB), intent(in) :: n
        integer (kind=IB) :: i,iter

        real    (kind=RP) :: r(2,n), x, f, df, dx

        real    (kind=RP), allocatable, dimension(:) :: p0, p1, tmp

        allocate(p0(1))
        allocate(p1(2))
        
        p0 = [1._RP]
        p1 = [1._RP, 0._RP]
        
        do k = 2,n
           tmp = ((2*k-1)*[p1,0._RP]-(k-1)*[0.0_RP, 0.0_RP,p0]) / k
           p0 = p1
           p1 = tmp
        end do

        do i = 1,n
           x = cos( pi*(i-0.25_RP)/(n+0.5_RP) )
           do iter = 1,10
              f = p1(1)
              df = 0.0_RP
              do k = 2,size(p1)
                 df = f + x*df
                 f  = p1(k) + x * f
              end do
              dx =  f / df
              x = x - dx
              if (abs(dx)<10*epsilon(dx)) exit
           end do
           r(1,i) = x
           r(2,i) = 2/((1-x**2)*df**2)
        end do

        deallocate(p0)
        deallocate(p1)
        
      END FUNCTION gaussquadrature
      
END SUBROUTINE GaussLegendreQuadrature
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_GaussLegendreQuadrature( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_GaussLegendreQuadrature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Fri Sep 30 15:30:49 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Number : 003                                                          !
!
! 1) GaussLegendreQuadrature within function gaussquadrature
! 2) author_GaussLegendreQuadrature



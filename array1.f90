program array1
    implicit none
    integer i
    integer, parameter :: nx = 20
    real(8) x(nx)
    real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
    real(8) :: dx=(xmax-xmin)/(nx-1)
    real(8) , dimension(nx):: a1,a2,a3,x1,b1
    x=[((i-1)*dx,i=1,nx)]
!-------------------------------------Tridiag section
    a1=1.0/dx**2
    a2=-2/dx**2
    a3=1.0/dx**2
    b1=sin(x)
    a1(1)=0.0
    a2(1)=1.0
    a3(1)=0
    a1(nx)=0.0
    a2(nx)=1.0
    a3(nx)=0
    print '(f10.5)',b1
    call solve_tridiag(a1,a2,a3,b1,x1,nx)
    print *, "---------------"
    print '(f10.5)',x1
end program array1
!--------------------------------------------------------------------------------------
! sub routine for triadiagonalsolver
!--------------------------------------------------------------------------------------
subroutine solve_tridiag(a,b,c,d,x,n)
    implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations

      integer,parameter :: r8 = kind(1.d0)

      integer,intent(in) :: n
      real(r8),dimension(n),intent(in) :: a,b,c,d
      real(r8),dimension(n),intent(out) :: x
      real(r8),dimension(n) :: cp,dp
      real(r8) :: m
      integer i

! initialize c-prime and d-prime
      cp(1) = c(1)/b(1)
      dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
       do i = 2,n
         m = b(i)-cp(i-1)*a(i)
         cp(i) = c(i)/m
         dp(i) = (d(i)-dp(i-1)*a(i))/m
       end do
! initialize x
       x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
      do i = n-1, 1, -1
        x(i) = dp(i)-cp(i)*x(i+1)
      end do

  end subroutine solve_tridiag

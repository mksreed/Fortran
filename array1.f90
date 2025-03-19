program array1
  implicit none
  real(16) err_max,err_avg
  integer i,nx
  nx=10
  do i=1,5
    call solve_1D_Poisson_2_Implicit(nx)
    call solve_1D_Poisson_2(nx)
    call solve_1D_Poisson_4(nx)
    nx=nx*2
  end do
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
!--------------------------------------------------------------------------------------
! subroutine for 1D Poisson 2nd order implicit
!--------------------------------------------------------------------------------------
subroutine solve_1D_Poisson_2_Implicit(nx)
  implicit none
  integer,intent(in) :: nx
  real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
  integer i
  real(8) x(nx)
  real(8) dx
  real(8) , dimension(nx) :: a1,a2,a3,x1,b1,ue,err
  real(16) err_max,err_avg
  dx=(xmax-xmin)/(nx-1)
  x=[((i-1)*dx,i=1,nx)]
!-------------------------------------call tridiag solver
  a1=1.0/dx**2
  a2=-2/dx**2
  a3=1.0/dx**2
  b1=sin(x)
  ue=-sin(x)
  a1(1)=0.0
  a2(1)=1.0
  a3(1)=0.0
  a1(nx)=0.0
  a2(nx)=1.0
  a3(nx)=0.0
! print '(f10.5)',b1
  call solve_tridiag(a1,a2,a3,b1,x1,nx)
  err=x1-ue
  err_avg=sum(err*err)/nx
  err_max=sqrt(maxval(err*err))
  !print *, "---------------"
  !  print '(2f10.5)',x1,err
  print '(3f16.11)',dx,err_max,err_avg
end subroutine solve_1D_Poisson_2_Implicit
!--------------------------------------------------------------------------------------
! subroutine for 1D Poisson 2nd order
!--------------------------------------------------------------------------------------
subroutine solve_1D_Poisson_2(nx)
  implicit none
  integer,intent(in) :: nx
  real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
  integer i
  real(8) x(nx)
  real(8) dx
  real(8) , dimension(nx) :: a1,a2,a3,x1,b1,ue,err
  real(16) err_max,err_avg
  dx=(xmax-xmin)/(nx-1)
  x=[((i-1)*dx,i=1,nx)]
!-------------------------------------c
  a1=1.0/dx**2
  a2=-2/dx**2
  a3=1.0/dx**2
  x1=0
  b1=sin(x)
  ue=-sin(x)
  do i=1,100000
    x1(2:nx-1)=(x1(3:nx)+x1(1:nx-2)-b1(2:nx-1)*dx**2)/2.
  end do
! print '(f10.5)',b1
  err=x1-ue
  err_avg=sum(err*err)/nx
  err_max=sqrt(maxval(err*err))
  !print *, "---------------"
  !  print '(2f10.5)',x1,err
  print '(3f16.11)',dx,err_max,err_avg
end subroutine solve_1D_Poisson_2
!--------------------------------------------------------------------------------------
! subroutine for 1D Poisson 4nd order
!--------------------------------------------------------------------------------------
subroutine solve_1D_Poisson_4(nx)
  implicit none
  integer,intent(in) :: nx
  real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
  integer i,it
  real(8) x(0:nx+1)
  real(8) dx
  real(8) , dimension(0:nx+1) :: x0,x1,b1,ue,err
  real(16) err_max,err_avg,err_max_prev
  dx=(xmax-xmin)/(nx-1)
  x=[((i-1)*dx,i=0,nx+1)]
!-------------------------------------call tridiag solver
  x1=0
  x0=x1
  b1=sin(x)
  ue=-sin(x)
  do it=1,100000
    !x0(2:nx-1)=(-x1(4:nx+1)+16*x1(3:nx)+16*x1(1:nx-2)-x1(0:nx-3)-12*b1(2:nx-1)*dx**2)/30.
    do i=2,nx-1
      x1(i)=(-x1(i+2)+16*x1(i+1)+16*x1(i-1)-x1(i-2)-b1(i)*12*dx*dx)/30.
    end do
    x1(0)=-x1(2)
    x1(nx+1)=-x1(nx-1)
    !x1=x0
    do i=0,nx+1
     !print '(2f16.11)',x1(i),ue(i)
    end do
    !print '(1f16.11)',it*1.0
  end do
  err=x1-ue
  err_avg=sum(err*err)/nx
  err_max=sqrt(maxval(err*err))
  print '(3f16.11)',dx,err_max,err_avg
end subroutine solve_1D_Poisson_4
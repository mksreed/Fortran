program periodic_tridiagonal_solver
    implicit none
    integer, parameter :: nx = 20 , n=5 ! Size of the matrix (change this value)
    double precision, dimension(n) :: a, b, c, d, x
    integer :: i
    
    ! Define the tridiagonal matrix coefficients and right-hand side vector
    data a /0.0, 1.0, 1.0, 1.0, 1.0/
    data b /5.0, 5.0, 5.0, 5.0, 5.0/
    data c /2.0, 2.0, 2.0, 2.0, 2.0/
    data d /1.0, 1.0, 1.0, 1.0, 1.0/

    double precision :: aa(nx),bb(nx),cc(nx),dd(nx),xx(nx),uu(nx)
    double precision :: dx
    dx=3.14159265358979323/(nx-1)
    xx=[((i-1)*dx,i=1,nx,1)]
    aa=1/dx**2
    cc=1/dx**2
    bb=-2/dx**2
    dd=sin(xx)
    uu=0.0
    !aa(nx)=cc(1)
    !cc(1)=aa(nx)
    !xx=[((i-1)*dx,i=1,nx)]

  
    ! Call the function to solve the system
    !call solve_periodic_tridiagonal(a, b, c, d, n, x)
    call solve_periodic_tridiagonal(aa, bb, cc, dd, nx, uu)
    
    ! Output the solution
    !print *, "Solution vector x:"
    !print *, x
    !write(*, '(F12.6)')x
    !write (*,'(4f12.6)')xx,uu,-sin(xx),uu+sin(xx)
    do i=1,nx
      write(*, '(3F12.6)')xx(i),uu(i),-sin(xx(i))
    end do
  
  contains
  
    subroutine solve_periodic_tridiagonal(a, b, c, d, n, x)
      double precision, dimension(n), intent(in) :: a, b, c, d
      integer, intent(in) :: n
      double precision, dimension(n), intent(out) :: x
      double precision, dimension(n) :: alpha, beta
      integer :: i
  
      ! Initialize alpha and beta
      alpha(1) = c(1) / b(1)
      beta(1) = d(1) / b(1)
      
      ! Forward sweep for Thomas algorithm with periodic boundary
      do i = 2, n
        alpha(i) = c(i) / (b(i) - a(i) * alpha(i-1))
        beta(i) = (d(i) - a(i) * beta(i-1)) / (b(i) - a(i) * alpha(i-1))
      end do
      
      ! Periodic correction for last row
      alpha(n) = c(n) / (b(n) - a(n) * alpha(n-1))
      beta(n) = (d(n) - a(n) * beta(n-1)) / (b(n) - a(n) * alpha(n-1))
      
      ! Back substitution with periodic boundary
      x(n) = (d(n) - a(n) * beta(n)) / (b(n) - a(n) * alpha(n-1))
      do i = n-1, 1, -1
        x(i) = beta(i) - alpha(i) * x(i+1)
      end do
      
    end subroutine solve_periodic_tridiagonal
  
  end program periodic_tridiagonal_solver

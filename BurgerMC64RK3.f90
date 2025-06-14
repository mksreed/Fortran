program viscous_burgers_rk3
  implicit none
  integer, parameter :: dp = kind(1.0d0), nx = 41
  real(dp), parameter :: L = 1.0d0, nu = 0.0001d0, dx = L / nx
  real(dp), parameter :: cfl = 0.4d0, tfinal = 20.8d0
  real(dp), dimension(-2:nx+3) :: u, u1, u2, rhs, u_exact, x
  integer :: i, n, nsteps
  real(dp) :: dt, dtv, dti, t = 0.0d0, umax, error

  ! Output file
    open(unit=10, file='outfile.txt')
    write(10,*) "X     U   UE  T"

  ! Initial condition: smooth sine wave
    nsteps=0
  do i = -2, nx+3
    x(i)=(i-1)*dx
    u(i) = sin(2.0d0 * 3.141592653589793d0 * (i-1)*dx)
    u(i)=2.*x(i)
    u_exact(i)=2.*x(i)
    write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
  end do

  ! Time stepping loop
  do while (t < tfinal)
    call apply_boundary_condition(u,nx)

    ! Estimate time step from convective CFL condition
    umax = maxval(abs(u(1:nx)))
    dti = cfl * dx / (umax + 1.0d-12)
    dtv = abs(umax)*dx**2/nu
    dt=min(dti,dtv)
    print *, dti,dtv,dt

    if (t + dt > tfinal) dt = tfinal - t

    ! ===== RK3 STAGE 1 =====
    rhs=0
    call compute_rhsb(u, rhs)
    u1 = u
    u1 = u + dt * rhs
    call apply_boundary_condition(u1,nx)

    ! ===== RK3 STAGE 2 =====
    call compute_rhsf(u1, rhs)
    !u2 = 0.75d0 * u + 0.25d0 * (u1 + dt * rhs)
    u2 = (u+u1)/2 + dt/2*rhs
    call apply_boundary_condition(u2,nx)
    u=u2

    ! ===== RK3 STAGE 3 =====
    !call compute_rhsb(u2, rhs)
    !u = (1.0d0/3.0d0) * u + (2.0d0/3.0d0) * (u2 + dt * rhs)
    !u=u1
    !call apply_boundary_condition(u)

    t = t + dt
    n = n + 1
    if (mod(n, 10) == 0) then
      print *, 't =', t, 'dt=', dt
      do i = 1, nx
        write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
      end do
    endif
  end do

  print *, 'Done. Final time:', t
  call exact_solution(u, u_exact, x, dx, nu, nx, t)
  !u_exact=u
    error=0
    do i = 1, nx
       write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
       error=error+(u(i)-u_exact(i))**2
    end do
    close(10)
    print *, dx, error

contains

  subroutine compute_rhs(u, rhs)
    real(dp), intent(in)  :: u(-2:nx+3)
    real(dp), intent(out) :: rhs(-2:nx+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine

    subroutine compute_rhsf(u, rhs)
    real(dp), intent(in)  :: u(-2:nx+3)
    real(dp), intent(out) :: rhs(-2:nx+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6f(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine

  subroutine compute_rhsb(u, rhs)
    real(dp), intent(in)  :: u(-2:nx+3)
    real(dp), intent(out) :: rhs(-2:nx+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6b(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine


  function dfdx6(f, i, dx) result(dfdx)
    real(dp), intent(in) :: f(-2:nx+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdx
    dfdx = (  -f(i-3) + 9*f(i-2) - 45*f(i-1) + 45*f(i+1) - 9*f(i+2) + f(i+3) ) / (60.0d0 * dx)
  end function

  function dfdx6f(f, i, dx) result(dfdxf)
    real(dp), intent(in) :: f(-2:nx+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdxf
    dfdxf = (  f(i+1)-f(i)) / (dx)
    dfdxf = (-f(i+2)+8.*f(i+1)-7.*f(i)) / (6.*dx)
    dfdxf = (f(i+3)-9.*f(i+2)+45.*f(i+1)-37.*f(i)) / (30.*dx)
  end function

  function dfdx6b(f, i, dx) result(dfdxb)
    real(dp), intent(in) :: f(-2:nx+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdxb
    dfdxb = (  f(i)-f(i-1)) / (dx)
    dfdxb = (f(i-2)-8.*f(i-1)+7.*f(i)) / (6.*dx)
    dfdxb = -(f(i-3)-9.*f(i-2)+45.*f(i-1)-37.*f(i)) / (30.*dx)
  end function

  function d2fdx2_4th(f, i, dx) result(d2fdx2)
    real(dp), intent(in) :: f(-2:nx+3), dx
    integer, intent(in)  :: i
    real(dp) :: d2fdx2
    d2fdx2 = ( -f(i-2) + 16*f(i-1) - 30*f(i) + 16*f(i+1) - f(i+2) ) / (12.0d0 * dx**2)
  end function

  subroutine apply_boundary_condition(f,nx)
    real(dp), intent(inout) :: f(-2:nx+3)
    integer :: nx
    !f(-2:0)     = f(nx-3:nx-1)
    !f(nx+1:n+3) = f(2:4)

      f(0)   =2.*f(1) -f(2)
      f(-1)  =2.*f(0) -f(1)
      f(-2)  =2.*f(-1) -f(0)
      f(nx+1)=2.*f(nx)-f(nx-1)
      f(nx+2)=2.*f(nx+1)-f(nx)
      f(nx+3)=2.*f(nx+2)-f(nx+1)
    
  end subroutine
  subroutine exact_solution(u, u_exact, x, dx, nu, nx, t)
      real(dp), intent(in) :: dx, nu, t
      real(dp), intent(inout) :: u(-2:nx+3), x(-2:nx+3)
      real(dp), intent(out) :: u_exact(-2:nx+3)
      integer :: i
      !real(dp), dimension(-2:nx+3) :: u_pred
      integer :: n, nx

      n = nx

     ! do i = 1, n-0
         u_exact=2*x/(1.+2*t)
      !end do
  end subroutine exact_solution
end program
program viscous_burgers_rk3
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: imax = 201
  real(dp), parameter :: L = 10.0d0, nu = 1.0001d0
  real(dp), parameter :: cfl = 0.4d0, tfinal = 1.1d0
  real(dp), dimension(-2:imax+3) :: u, u1, u2, rhs, u_exact, x
  integer :: i,  nerror=0, nx, nsteps, file_unit=10
  real(dp) :: dt, dtv, dti, t = 0.0d0, umax, error(10),dx, error_slope
  real(dp) :: dxs(10)

  ! Output file
    open(unit=10, file='outfile.txt')
    write(file_unit,*) "X     U   UE  T"

    do nx= 41,161,40
  ! Initial condition: smooth sine wave
    nsteps=0
    !nx=41
    dx=2.*L/nx
    u=0
    u_exact=0
    x=0
    t=1
    call initialize(u,u_exact,x,nx,imax, dx, nu, t)
    do i = 1, nx
     write(file_unit,'(4E15.7)') x(i), u(i), u_exact(i),t
    end do

  ! Time stepping loop
    call timeintegration(u,x, dx,cfl,t,tfinal,nx,imax)
  

  !print *, 'Done. Final time:', t
  call exact_solution(u, u_exact, x, dx, nu, nx,imax, t)
  !u_exact=u
    nerror=nerror+1
    error(nerror)=0
    dxs(nerror)=dx
    do i = 1, nx
       write(file_unit,'(4E15.7)') x(i), u(i), u_exact(i),t
       error(nerror)=error(nerror)+(u(i)-u_exact(i))**2
    end do

    print *, nx, dx, error(nerror),t
  end do
  do i=2,nerror
    error_slope=(log(error(i))-log(error(i-1)))/(log(dxs(i))-log(dxs(i-1)))
    print *, i, dxs(i),error(i),error_slope
  end do

  close(10)

contains
 subroutine timeintegration(u,x,dx,cfl,t,tfinal,nx,imax)
  integer :: nx, imax
  real(dp)  ::  dx
  real(dp) :: cfl , tfinal
  real(dp), dimension(-2:imax+3) :: u, u1, u2, rhs, u_exact, x
  integer :: i, n, nsteps
  real(dp) :: dt, dtv, dti, t
  do while (t < tfinal)
    call apply_boundary_condition(u,nx,imax)

    ! Estimate time step from convective CFL condition
    umax = maxval(abs(u(1:nx)))
    dti = cfl * dx / (umax + 1.0d-12)
    dtv = abs(umax)*dx**2/nu
    dt=min(dti,dtv)
   ! print *, n, dti,dtv,dt

    if (t + dt > tfinal) dt = tfinal - t

    ! ===== RK3 STAGE 1 =====
    rhs=0
    call compute_rhsb(u, rhs)
    u1 = u
    u1 = u + dt * rhs
    call apply_boundary_condition(u1,nx,imax)

    ! ===== RK3 STAGE 2 =====
    call compute_rhsf(u1, rhs)
    !u2 = 0.75d0 * u + 0.25d0 * (u1 + dt * rhs)
    u2 = (u+u1)/2 + dt/2*rhs
    call apply_boundary_condition(u2,nx,imax)
    u=u2

    ! ===== RK3 STAGE 3 =====
    !call compute_rhsb(u2, rhs)
    !u = (1.0d0/3.0d0) * u + (2.0d0/3.0d0) * (u2 + dt * rhs)
    !u=u1
    !call apply_boundary_condition(u)

    t = t + dt
    n = n + 1
    if (mod(n, 5000000) == 0) then
      print *, 't =', t, 'dt=', dt
      do i = 1, nx
        write(file_unit,'(4E15.7)') x(i), u(i), u_exact(i),t
      end do
    endif
  end do

 end subroutine timeintegration

  subroutine initialize(u,u_exact,x,nx,imax,dx,nu,t)
    integer :: i, nx, imax
    real(dp), intent(inout)  :: u(-2:imax+3), x(-2:imax+3), u_exact(-2:imax+3)
    real (dp), intent(in) :: dx, t, nu
    real(dp) :: pi
  ! Initial condition: smooth sine wave
    pi=4.*atan(1.0)
    nsteps=0
   do i = -2, nx+3
    x(i)=(i-nx/2)*dx
    u(i) = sin(2.0d0 * 3.141592653589793d0 * (i-1)*dx)
    u(i)=2.*x(i)
    u_exact(i)=2.*x(i)
   end do
   u=x/t-pi/t*tanh(pi*x/2./nu/t)
   u_exact=x/t-pi/t*tanh(pi*x/2./nu/t)
  end subroutine initialize

  subroutine compute_rhs(u, rhs)
    real(dp), intent(in)  :: u(-2:imax+3)
    real(dp), intent(out) :: rhs(-2:imax+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine

    subroutine compute_rhsf(u, rhs)
    real(dp), intent(in)  :: u(-2:imax+3)
    real(dp), intent(out) :: rhs(-2:imax+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6f(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine

  subroutine compute_rhsb(u, rhs)
    real(dp), intent(in)  :: u(-2:imax+3)
    real(dp), intent(out) :: rhs(-2:imax+3)
    integer :: i
    real(dp) :: dudx, d2udx2

    do i = 1, nx
       dudx    = dfdx6b(u, i, dx)
       d2udx2  = d2fdx2_4th(u, i, dx)
       rhs(i)  = -u(i) * dudx + nu * d2udx2
    end do
  end subroutine


  function dfdx6(f, i, dx) result(dfdx)
    real(dp), intent(in) :: f(-2:imax+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdx
    dfdx = (  -f(i-3) + 9*f(i-2) - 45*f(i-1) + 45*f(i+1) - 9*f(i+2) + f(i+3) ) / (60.0d0 * dx)
  end function

  function dfdx6f(f, i, dx) result(dfdxf)
    real(dp), intent(in) :: f(-2:imax+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdxf
    dfdxf = (  f(i+1)-f(i)) / (dx)
    dfdxf = (-f(i+2)+8.*f(i+1)-7.*f(i)) / (6.*dx)
    dfdxf = (f(i+3)-9.*f(i+2)+45.*f(i+1)-37.*f(i)) / (30.*dx)
  end function

  function dfdx6b(f, i, dx) result(dfdxb)
    real(dp), intent(in) :: f(-2:imax+3), dx
    integer, intent(in)  :: i
    real(dp) :: dfdxb
    dfdxb = (  f(i)-f(i-1)) / (dx)
    dfdxb = (f(i-2)-8.*f(i-1)+7.*f(i)) / (6.*dx)
    dfdxb = -(f(i-3)-9.*f(i-2)+45.*f(i-1)-37.*f(i)) / (30.*dx)
  end function

  function d2fdx2_4th(f, i, dx) result(d2fdx2)
    real(dp), intent(in) :: f(-2:imax+3), dx
    integer, intent(in)  :: i
    real(dp) :: d2fdx2
    d2fdx2 = ( -f(i-2) + 16*f(i-1) - 30*f(i) + 16*f(i+1) - f(i+2) ) / (12.0d0 * dx**2)
  end function

  subroutine apply_boundary_condition(f,nx, imax)
    integer :: imax
    real(dp), intent(inout) :: f(-2:imax+3)
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
  subroutine exact_solution(u, u_exact, x, dx, nu, nx, imax, t)
      real(dp), intent(in) :: dx, nu, t
      integer :: imax
      real(dp), intent(inout) :: u(-2:imax+3), x(-2:imax+3)
      real(dp), intent(out) :: u_exact(-2:imax+3)
      integer :: i
      real (dp) :: pi
      !real(dp), dimension(-2:imax+3) :: u_pred
      integer :: n, nx
      pi=4.*atan(1.0)

      n = nx

     ! do i = 1, n-0
         u_exact=2*x/(1.+2*t)
         u_exact=x/t-pi/t*tanh(pi*x/2./nu/t)
      !end do
  end subroutine exact_solution
end program
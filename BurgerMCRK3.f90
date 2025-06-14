program burgers_rk3_maccormack
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nx = 101
    real(dp), parameter :: L = 1.0_dp
    real(dp), parameter :: nu = 0.0001_dp   ! viscosity
    real(dp), parameter :: CFL = 0.3_dp
    real(dp), parameter :: t_end = 1.8_dp
  
    real(dp) :: dx, dt, t, error
    real(dp), dimension(-2:nx+3) :: x, u, u1, u2, rhs, u_exact
  
    integer :: i, nsteps

    open(unit=10, file='outfile.txt')
    write(10,*) "X     U   UE  T"
    dx = L / (nx - 1)
    do i = 0, nx+1
       x(i) = (i - 1) * dx
    end do
  
    ! Initial condition: smooth step
    u = 1.0_dp
    do i = 0, nx+1
       if (x(i) > 0.5_dp) u(i) = 0.0_dp
       u(i)=abs(sin(4*3.14159265358797323*x(i)))
       u(i)=2.*x(i)
    end do
  
    t = 0.0_dp
    nsteps = 0
    dt = CFL * dx / maxval(abs(u))
    call exact_solution(u, u_exact, x, dx, nu, nx, t)
   do i = 1, nx
      write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
   end do
  
    do while (t < t_end)
       if (t + dt > t_end) dt = t_end - t
  
       ! RK stage 1
       call compute_rhs(u, rhs, dx, nu, nx)
       u1 = u - dt * rhs
  
       ! RK stage 2
       call compute_rhs(u1, rhs, dx, nu, nx)
       u2 = 0.75_dp * u + 0.25_dp * (u1 - dt * rhs)
  
       ! RK stage 3
       call compute_rhs(u2, rhs, dx, nu, nx)
       u = (1.0_dp/3.0_dp) * u + (2.0_dp/3.0_dp) * (u2 - dt * rhs)
  
       t = t + dt
       nsteps = nsteps + 1

       if(mod(nsteps,100)==0) then
         call exact_solution(u, u_exact, x, dx, nu, nx, t)
          do i = 1, nx
             write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
          end do
       end if
    end do
  
    ! Output results
    call exact_solution(u, u_exact, x, dx, nu, nx, t)
    error=0
    do i = 1, nx
       write(10,'(4E15.7)') x(i), u(i), u_exact(i),t
       error=error+(u(i)-u_exact(i))**2
    end do
    close(10)
    print *, dx, error
    print *, 'Simulation complete. Steps:', nsteps, ' Final time:', t
  
  contains
  
    subroutine compute_rhs(u, rhs, dx, nu, nx)
      real(dp), intent(in) :: dx, nu
      real(dp), intent(inout) :: u(-2:nx+3)
      real(dp), intent(out) :: rhs(-2:nx+3)
      integer :: i
      real(dp), dimension(-2:nx+3) :: u_pred, dudx, d2udx2, f1
      integer :: n, nx
  
      n = size(u)
      n = nx

      ! Apply Periodic BC
      u(-2:0)=u(nx-3:nx-1)
      u(n+1:n+3)=u(2:4)

      !Apply Linear extrapolation
      u(0)   =2.*u(1) -u(2)
      u(-1)  =2.*u(0) -u(1)
      u(nx+1)=2.*u(nx)-u(nx-1)
      u(nx+2)=2.*u(nx+1)-u(nx)

      !computing flux
      f1=0.5*u(1)**2

      ! Predictor step (forward)

      u_pred = u
      do i = 1, n-0
         u_pred(i) = u(i) - dx * 0.5_dp * u(i) * (u(i+1) - u(i)) / dx
         !u_pred(i) = u(i) - dt * (-f1(i+2)+8.*f1(i+1) - 7.*f1(i)) / (6*dx)
      end do
  
      ! Corrector step (backward)
      dudx = 0.0_dp
      do i = 1, n-0
         dudx(i) = 0.5_dp * ((u_pred(i) * (u_pred(i) - u_pred(i-1)) / dx) + &
                             (u(i) * (u(i+1) - u(i)) / dx))
      end do
  
      ! Second-order central diff for viscous term
      d2udx2 = 0.0_dp
      do i = 1, n-0
         d2udx2(i) = (u(i+1) - 2.0_dp*u(i) + u(i-1)) / dx**2
      end do
  
      rhs = dudx - nu * d2udx2
  
      ! Boundary conditions: Dirichlet u = 1 (left), u = 0 (right)
      !rhs(1) = 0.0_dp
      !rhs(n) = 0.0_dp
      !rhs(0) = 0.0_dp
      !rhs(n+1) = 0.0_dp
    end subroutine
    !--------------------------------------------------
    subroutine exact_solution(u, u_exact, x, dx, nu, nx, t)
      real(dp), intent(in) :: dx, nu, t
      real(dp), intent(inout) :: u(-2:nx+3), x(-2:nx+3)
      real(dp), intent(out) :: u_exact(-2:nx+3)
      integer :: i
      !real(dp), dimension(-2:nx+3) :: u_pred
      integer :: n, nx
  
      n = size(u)
      n = nx

      do i = 1, n-0
         u_exact=2*x/(1.+2*t)
      end do
  end subroutine exact_solution
  end program burgers_rk3_maccormack

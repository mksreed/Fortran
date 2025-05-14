program burgers_rk3_maccormack
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nx = 101
    real(dp), parameter :: L = 1.0_dp
    real(dp), parameter :: nu = 0.01_dp   ! viscosity
    real(dp), parameter :: CFL = 0.4_dp
    real(dp), parameter :: t_end = 0.5_dp
  
    real(dp) :: dx, dt, t
    real(dp), dimension(nx) :: x, u, u1, u2, rhs
  
    integer :: i, nsteps

    open(unit=10, file='outfile.txt')
    write(10,*) "X     U    T"
    dx = L / (nx - 1)
    do i = 1, nx
       x(i) = (i - 1) * dx
    end do
  
    ! Initial condition: smooth step
    u = 1.0_dp
    do i = 1, nx
       if (x(i) > 0.5_dp) u(i) = 0.0_dp
    end do
  
    t = 0.0_dp
    nsteps = 0
    dt = CFL * dx / maxval(abs(u))
  
    do while (t < t_end)
       if (t + dt > t_end) dt = t_end - t
  
       ! RK stage 1
       call compute_rhs(u, rhs, dx, nu)
       u1 = u - dt * rhs
  
       ! RK stage 2
       call compute_rhs(u1, rhs, dx, nu)
       u2 = 0.75_dp * u + 0.25_dp * (u1 - dt * rhs)
  
       ! RK stage 3
       call compute_rhs(u2, rhs, dx, nu)
       u = (1.0_dp/3.0_dp) * u + (2.0_dp/3.0_dp) * (u2 - dt * rhs)
  
       t = t + dt
       nsteps = nsteps + 1

       if(mod(nsteps,20)==0) then
          do i = 1, nx
             write(10,'(3E15.7)') x(i), u(i), t
          end do
       end if
    end do
  
    ! Output results

    do i = 1, nx
       write(10,'(3E15.7)') x(i), u(i), t
    end do
    close(10)
  
    print *, 'Simulation complete. Steps:', nsteps, ' Final time:', t
  
  contains
  
    subroutine compute_rhs(u, rhs, dx, nu)
      real(dp), intent(in) :: u(:), dx, nu
      real(dp), intent(out) :: rhs(size(u))
      integer :: i
      real(dp), dimension(size(u)) :: u_pred, dudx, d2udx2
      integer :: n
  
      n = size(u)
  
      ! Predictor step (forward)
      u_pred = u
      do i = 2, n-1
         u_pred(i) = u(i) - dx * 0.5_dp * u(i) * (u(i+1) - u(i)) / dx
      end do
  
      ! Corrector step (backward)
      dudx = 0.0_dp
      do i = 2, n-1
         dudx(i) = 0.5_dp * ((u_pred(i) * (u_pred(i) - u_pred(i-1)) / dx) + &
                             (u(i) * (u(i+1) - u(i)) / dx))
      end do
  
      ! Second-order central diff for viscous term
      d2udx2 = 0.0_dp
      do i = 2, n-1
         d2udx2(i) = (u(i+1) - 2.0_dp*u(i) + u(i-1)) / dx**2
      end do
  
      rhs = dudx - nu * d2udx2
  
      ! Boundary conditions: Dirichlet u = 1 (left), u = 0 (right)
      rhs(1) = 0.0_dp
      rhs(n) = 0.0_dp
    end subroutine
  
  end program burgers_rk3_maccormack
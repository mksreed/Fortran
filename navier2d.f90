program navier_stokes_2d
    implicit none
    integer, parameter :: nx = 50, ny = 50, nt = 500
    double precision, parameter :: dx = 1.0/nx, dy = 1.0/ny
    double precision, parameter :: dt = 0.001, Re = 100.0
    double precision, dimension(0:nx,0:ny) :: x, y, u, v, u_new, v_new, p,  b
    double precision :: umax=0,vmax=0,theta=0,rad=0
    integer :: i, j, n
    double precision :: rho = 1.0, nu, time=0
    character(len=20) :: filename, filename1
    integer :: unit=10, unit1=11
    filename = "outfileuv.txt"
    filename1 = "outfile.txt"
    open(unit=unit, file=filename, status='replace')
    open(unit=unit1, file=filename1, status='replace')
    write(unit1, '(A)') "t,umax,vmax"
    write(unit, '(A)') "x,y,u,v,t"

    nu = 1.0 / Re
  
    ! Initialize fields
    u = 1.0
    v = 0.0
    p = 0.0
    do j=0, ny
      do i = 0,nx
         x(i,j)=(i-nx/2)*dx
         y(i,j)=(j-ny/2)*dy
         theta=atan(y(i,j),x(i,j))
         rad=sqrt(x(i,j)**2+y(i,j)**2)
         u(i,j)=0.-rad*sin(theta)*exp(-5.*rad)*2
         v(i,j)=rad*cos(theta)*exp(-5.*rad)*2
      end do
     end do
  
    do n = 1, nt
      call build_rhs(b, u, v, dx, dy, dt, rho,nx,ny)
      call pressure_poisson(p, b, dx, dy,nx ,ny)
      call update_velocity(u, v, p, u_new, v_new, dx, dy, dt, rho, nu,nx ,ny)
      umax=maxval(u-u_new)
      vmax=maxval(v-v_new)
      u = u_new
      v = v_new
      time=time+dt

      write(unit1,'(3F12.6)') time, umax, vmax
      write(unit1,'(10F12.6)') u_new(15,:)
    end do
  
    print *, 'Simulation finished.'
    do j=0, ny
     do i = 0,nx
        write(unit,'(5F12.6)') x(i,j),y(i,j),u(i,j),v(i,j),time
     end do
    end do

  end program navier_stokes_2d
  !----------------------------------------------
  subroutine build_rhs(b, u, v, dx, dy, dt, rho,nx,ny)
    implicit none
    double precision, intent(in) :: dx, dy, dt, rho
    double precision, dimension(0:nx,0:ny), intent(in) :: u, v
    double precision, dimension(0:nx,0:ny), intent(out) :: b
    integer :: i, j, nx, ny
    !nx = size(u,1)-1
    !ny = size(u,2)-1
  
    do i = 1, nx-1
      do j = 1, ny-1
        b(i,j) = rho * (1.0/dt * ((u(i+1,j) - u(i-1,j))/(2*dx) + (v(i,j+1) - v(i,j-1))/(2*dy)))
      end do
    end do
  end subroutine build_rhs
!-------------------------------------------------
  subroutine pressure_poisson(p, b, dx, dy,nx,ny)
    implicit none
    double precision, intent(in) :: dx, dy
    double precision, dimension(0:nx,0:ny), intent(in) :: b
    double precision, dimension(0:nx,0:ny), intent(inout) :: p
    double precision :: pn(0:size(p,1)-1,0:size(p,2)-1)
    integer :: i, j, n, nx, ny
    !nx = size(p,1)-1
    !ny = size(p,2)-1
  
    do n = 1, 50
      pn = p
      do i = 1, nx-1
        do j = 1, ny-1
          p(i,j) = (((pn(i+1,j) + pn(i-1,j)) * dy**2 + &
          (pn(i,j+1) + pn(i,j-1)) * dx**2) - b(i,j) * dx**2 * dy**2) &
          / (2.0 * (dx**2 + dy**2))
        end do
      end do
  
      ! Boundary conditions: Neumann (dp/dn = 0)
      p(:,ny) = p(:,ny-1)
      p(:,0)  = p(:,1)
      p(0,:)  = p(1,:)
      p(nx,:) = p(nx-1,:)
    end do
  end subroutine pressure_poisson
!----------------------------------------------------------
  subroutine update_velocity(u, v, p, u_new, v_new, dx, dy, dt, rho, nu,nx,ny)
    implicit none
    double precision, intent(in) :: dx, dy, dt, rho, nu
    double precision, dimension(0:nx,0:ny), intent(in) :: u, v, p
    double precision, dimension(0:nx,0:ny), intent(out) :: u_new, v_new
    integer :: i, j, nx, ny
  
    !nx = size(u,1)-1
    !ny = size(u,2)-1
  
    do i = 1, nx-1
      do j = 1, ny-1
        u_new(i,j) = u(i,j) - u(i,j)*(dt/dx)*(u(i,j) - u(i-1,j)) - v(i,j)*(dt/dy)*(u(i,j) - u(i,j-1)) &
                    - dt/(2*rho*dx)*(p(i+1,j) - p(i-1,j)) + nu*(dt/dx**2)*(u(i+1,j)-2*u(i,j)+u(i-1,j)) &
                    + nu*(dt/dy**2)*(u(i,j+1)-2*u(i,j)+u(i,j-1))
  
        v_new(i,j) = v(i,j) - u(i,j)*(dt/dx)*(v(i,j) - v(i-1,j)) - v(i,j)*(dt/dy)*(v(i,j) - v(i,j-1)) &
                    - dt/(2*rho*dy)*(p(i,j+1) - p(i,j-1)) + nu*(dt/dx**2)*(v(i+1,j)-2*v(i,j)+v(i-1,j)) &
                    + nu*(dt/dy**2)*(v(i,j+1)-2*v(i,j)+v(i,j-1))
      end do
    end do
  
    ! Boundary Conditions (no-slip) Drv Cvt
    u_new(0,:) = 0.0; u_new(nx,:) = 0.0
    u_new(:,0) = 0.0; u_new(:,ny) = 1.0
    v_new(0,:) = 0.0; v_new(nx,:) = 0.0
    v_new(:,0) = 0.0; v_new(:,ny) = 0.0
        ! Boundary Conditions (no-slip) BL
    u_new(0,:) = 1.0*0.0; u_new(nx,:) = u_new(nx-1,:)*0
    u_new(:,0) = 0.0; u_new(:,ny) = 1.0*0
    v_new(0,:) = 0.0; v_new(nx,:) = v_new(nx-1,:)*0
    v_new(:,0) = 0.0; v_new(:,ny) = 0.0
  end subroutine update_velocity
      
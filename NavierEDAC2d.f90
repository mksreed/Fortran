program edac_navier_stokes
    implicit none
    integer, parameter :: nx = 25, ny = 25, nt = 1000
    double precision, parameter  :: pi=3.14159265358979323
    double precision, parameter :: dx = 2.*pi/(nx-1), dy = 2.*pi/(ny-1), dt = 0.001
    double precision, parameter :: nu = 0.01, c = 1.0 
    double precision time
  
    double precision :: u(nx,ny), v(nx,ny), p(nx,ny),x(nx,ny),y(nx,ny)
    double precision :: un(nx,ny), vn(nx,ny), pn(nx,ny)
    double precision :: ue(nx,ny), ve(nx,ny), pe(nx,ny)
    double precision :: umax,vmax,pmax
    integer :: i, j, n
    double precision :: dudx, dudy, dvdx, dvdy, div, lap_p
    character(len=20) :: filename, filename1

    integer :: unit=10, unit1=11
    filename = "outfileuv.txt"
    filename1 = "outfile.txt"
    open(unit=unit, file=filename, status='replace')
    open(unit=unit1, file=filename1, status='replace')
    write(unit1, '(A)') "t,umax,vmax,pmax"
    write(unit, '(A)') "x,y,u,v,t"
  
    ! Initialize fields
    u = 0.0; v = 0.0; p = 0.0; time=0
    do j=1, ny
        do i = 1,nx
           x(i,j)=(i-nx/2*0)*dx
           y(i,j)=(j-ny/2*0)*dy
           !theta=atan(y(i,j),x(i,j))
           !rad=sqrt(x(i,j)**2+y(i,j)**2)
           !u(i,j)=0.-rad*sin(theta)*exp(-5.*rad)*2
           !v(i,j)=rad*cos(theta)*exp(-5.*rad)*2
           x(i,j)=(i-1)*dx
           y(i,j)=(j-1)*dy
        end do
       end do
       call taylor_vortices(ue,ve,pe,x,y,nx,ny,nu,time)
       u=ue
       v=ve
       p=pe
  
    do n = 1, nt
       un = u
       vn = v
       pn = p
  
       ! Update velocity (u, v)
       do i = 2, nx-1
          do j = 2, ny-1
             dudx = (pn(i+1,j) - pn(i-1,j)) / (2.0*dx)
             dudy = (un(i,j+1) - un(i,j-1)) / (2.0*dy)
             dvdx = (vn(i+1,j) - vn(i-1,j)) / (2.0*dx)
             dvdy = (pn(i,j+1) - pn(i,j-1)) / (2.0*dy)
  
             u(i,j) = un(i,j) + dt * ( &
       &         - un(i,j) * (un(i+1,j) - un(i-1,j)) / (2.0*dx) &
       &         - vn(i,j) * (un(i,j+1) - un(i,j-1)) / (2.0*dy) &
       &         - dudx &
       &         + nu * ((un(i+1,j) - 2.0*un(i,j) + un(i-1,j)) / dx**2 + &
       &                 (un(i,j+1) - 2.0*un(i,j) + un(i,j-1)) / dy**2) )
  
             v(i,j) = vn(i,j) + dt * ( &
       &         - un(i,j) * (vn(i+1,j) - vn(i-1,j)) / (2.0*dx) &
       &         - vn(i,j) * (vn(i,j+1) - vn(i,j-1)) / (2.0*dy) &
       &         - dvdy &
       &         + nu * ((vn(i+1,j) - 2.0*vn(i,j) + vn(i-1,j)) / dx**2 + &
       &                 (vn(i,j+1) - 2.0*vn(i,j) + vn(i,j-1)) / dy**2) )
          end do
       end do
  
       ! EDAC pressure update
       do i = 2, nx-1
          do j = 2, ny-1
             div = ((u(i+1,j) - u(i-1,j)) / (2.0*dx) + &
       &            (v(i,j+1) - v(i,j-1)) / (2.0*dy)) 
  
             lap_p = (pn(i+1,j) - 2.0*pn(i,j) + pn(i-1,j)) / dx**2 + &
       &             (pn(i,j+1) - 2.0*pn(i,j) + pn(i,j-1)) / dy**2
  
             p(i,j) = pn(i,j) - dt * (c**2 * div - nu * lap_p)
          end do
       end do

       time=time+dt
       call taylor_vortices(ue,ve,pe,x,y,nx,ny,nu,time)
       ! Boundary conditions
       u(:,ny) = 1.0     ! Lid driven at top
       u(:,1)  = 0.0; u(1,:) = 0.0; u(nx,:) = 0.0
       v(:,ny) = 0.0
       v(:,1)  = 0.0; v(1,:) = 0.0; v(nx,:) = 0.0
       p(1,:) = p(2,:); p(nx,:) = p(nx-1,:)
       p(:,1) = p(:,2); p(:,ny) = p(:,ny-1)
       !taylor vortex
       u(:,ny) = ue(:,ny)     
       u(:,1)  = ue(:,1); u(1,:) = ue(1,:); u(nx,:) = ue(nx,:) 
       v(:,ny) = ve(:,ny)
       v(:,1)  = ve(:,1) ; v(1,:) = ve(1,:); v(nx,:) = ve(nx,:) 
       p(1,:) = pe(1,:); p(nx,:) = pe(nx,:)
       p(:,1) = pe(:,1); p(:,ny) = pe(:,ny) 
       !
       umax=maxval(u-ue)
       vmax=maxval(v-ve)
       pmax=maxval(p-pe)
       write(unit1,'(4F12.6)') time, umax,vmax,pmax
    end do


    print *, "Simulation completed."
    do j=1, ny
        do i = 1,nx
           write(unit,'(5F12.6)') x(i,j),y(i,j),u(i,j),v(i,j),time
        end do
       end do
  end program edac_navier_stokes
  !***********************************************************************
  subroutine taylor_vortices(u,v,p,x,y,nx,ny,nu,t)
    integer, intent(in)  :: nx, ny 
    double precision :: u(nx,ny), v(nx,ny), p(nx,ny),x(nx,ny),y(nx,ny)
    double precision nu,t
    do j = 1, ny
        do i = 1, nx
           u(i,j) = -cos(x(i,j)) * sin(y(i,j)) * exp(-2.0 * nu * t)
           v(i,j) =  sin(x(i,j)) * cos(y(i,j)) * exp(-2.0 * nu * t)
           p(i,j) = -0.25 * (cos(2.0 * x(i,j)) + cos(2.0 * y(i,j))) * exp(-4.0 * nu * t)
        end do
     end do
  end subroutine  taylor_vortices
  
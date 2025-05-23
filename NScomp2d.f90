program navier_stokes_fd
    implicit none
    ! Parameters
    integer, parameter :: nx = 40, ny = 60
    integer, parameter :: imax = nx+3, jmax = ny+3
    integer, parameter :: nt = 10000
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax=1, ymax=0.2
    double precision, parameter :: dt = 1e-6
    double precision, parameter :: gamma = 1.4, mach=0.6, s1=110.4/273., re=100, pr=0.72
    double precision, parameter :: R = 287.0
    double precision, parameter :: mu = 1.8e-5
    double precision, parameter :: kappa = 0.025
    integer i,j


    ! State variables
    double precision :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
    double precision :: p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
    double precision :: rho_new(-2:imax,-2:jmax), u_new(-2:imax,-2:jmax), v_new(-2:imax,-2:jmax)
    double precision ::  T_new(-2:imax,-2:jmax)
    double precision :: time

    integer :: n

    ! Initialization
    call initialize(x,y,dx,dy,xmax,ymax,rho, u, v, T,nx,ny,imax,jmax)
    call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
    call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny,imax,jmax)
    !do j=1,ny,1; do i=1,nx,1;write(*,'(2I4,3F20.13)'),i,j,u(i,j),v(i,j),p(i,j);end do;end do
    rho_new = rho
    u_new = u
    v_new = v
    T_new = T
    ! Time-stepping
    time=0
      do n = 1, nt
        call RK2_update(rho, u, v, T, p, E,rho_new,u_new,v_new,T_new,re,pr,gamma,mach, &
                                                    s1,dt,dx,dy,nx,ny,imax,jmax)
        call apply_bc(rho_new, u_new, v_new, T_new,nx,ny,imax,jmax)
        !call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
        call compute_pressure_energy(rho_new,u_new,v_new,T_new,p,E,re,pr,gamma,mach, &
                                                                s1,nx,ny,imax,jmax)

        rho = rho_new
        u = u_new
        v = v_new
        T = T_new
        time=time+dt
        write(*,'(6F12.6)'),time,maxval(rho),maxval(u),maxval(v),maxval(T),maxval(p)
        write(*,'(4I5)') maxloc(rho),maxloc(p)
    end do

    call write_output(x,y,rho, u, v, T, p,nx,ny,imax,jmax,time)
end program navier_stokes_fd
!---------------------
subroutine initialize(x,y,dx,dy,xmax,ymax,rho, u, v, T,nx,ny,imax,jmax)
    implicit none
    double precision, intent(out) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), v(-2:imax,-2:jmax), &
                                                 T(-2:imax,-2:jmax)
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax, ymax
    double precision :: rad, theta, pi


    integer :: i, j, nx, ny, imax, jmax

    dx=xmax/(nx-1)
    dy=ymax/(ny-1)
    pi=4.*atan(45.0)
    do i=-2,imax,1;x(i)=(i-1)*dx;end do
    do j=-2,jmax,1;y(j)=(j-1)*dy; end do

    do j = -2, jmax
        do i = -2, imax
            rho(i,j) = 1.0
            u(i,j) = tanh(40*y(j))+sin(200*x(i)*pi)/20.0*exp(-y(j)*y(j)*1000)
            u(i,j)=1.0
            v(i,j) = 0.0
            T(i,j) = 1.0
            theta=atan(y(j),x(i))
            rad=sqrt((x(i))**2+(y(j))**2)
            !u(i,j)=u(i,j)+0.-rad*sin(theta)*exp(-5.*rad)*10
            !v(i,j)=v(i,j)+rad*cos(theta)*exp(-5.*rad)*10
            !write(*,'(2I4,4F20.13)'),i,j,u(i,j),v(i,j),rad,theta
        end do
    end do
end subroutine initialize
!--------------------------------------
subroutine apply_bc(rho, u, v, T,nx,ny,imax,jmax)
    implicit none
    double precision, intent(inout) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
                                        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
    double precision :: r1(nx),u1(nx),v1(nx),t1(nx)
    double precision :: rf(nx),uf(nx),vf(nx),tf(nx)
    double precision :: dx,dy
    integer :: i, j, nx, ny, imax, jmax
    integer :: bp=0,b0=1,bs=0,bn=0,ba=0;
    double precision :: temp

    do j = 1, ny
        rho(1,j) = rho(2,j)*b0+rho(nx,j)*bp
        u(1,j)   =   u(2,j)*b0+  u(nx,j)*bp
        v(1,j)   =   v(2,j)*b0+  v(nx,j)*bp
        T(1,j)   =   T(2,j)*b0+  T(nx,j)*bp
        rho(0,j) = rho(2,j)*b0+rho(nx-1,j)*bp
        u(0,j)   =   u(2,j)*b0+  u(nx-1,j)*bp
        v(0,j)   =   v(2,j)*b0+  v(nx-1,j)*bp
        T(0,j)   =   T(2,j)*b0+  T(nx-1,j)*bp
        u(-2:1,j)=1.0
        v(-2:1,j)=0.
        T(-2:1,j)=1.0
        rho(-2:1,j)=1.0

        rho(nx,j) = rho(nx-1,j)*b0+rho(1,j)*bp
        u(nx,j)   =   u(nx-1,j)*b0+  u(1,j)*bp
        v(nx,j)   =   v(nx-1,j)*b0+  v(1,j)*bp
        T(nx,j)   =   T(nx-1,j)*b0+  T(1,j)*bp
        rho(nx+1,j) = rho(nx-1,j)*b0+rho(2,j)*bp
        u(nx+1,j)   =   u(nx-1,j)*b0+  u(2,j)*bp
        v(nx+1,j)   =   v(nx-1,j)*b0+  v(2,j)*bp
        T(nx+1,j)   =   T(nx-1,j)*b0+  T(2,j)*bp
        !u(nx,j)=0.
        !v(nx,j)=0.
        !u(nx+1,j)=0.
        !v(nx+1,j)=0.
    end do

    do i = 1, nx
        rho(i,1) = rho(i,2)*b0+rho(i,ny)*bp
        u(i,1)   =   u(i,2)*b0+  u(i,ny)*bp
        v(i,1)   =   v(i,2)*b0+  v(i,ny)*bp
        T(i,1)   =   T(i,2)*b0+  T(i,ny)*bp
        rho(i,0) = rho(i,2)*b0+rho(i,ny-1)*bp
        u(i,0)   =   u(i,2)*b0+  u(i,ny-1)*bp
        v(i,0)   =   v(i,2)*b0+  v(i,ny-1)*bp
        T(i,0)   =   T(i,2)*b0+  T(i,ny-1)*bp
        u(i,-2:1)=0  
        v(i,-2:1)=0
        T(i,-2:1)=T(i,2)
        rho(i,-2:1)=rho(i,2)

        !do i = 1, nx-0,5
        !do j = 1, ny-0
            !write(12,'(2I4,4E14.7)'),i,j, fr(i,j),fu(i,j),fv(i,j),fe(i,j)!madhu
            !write(12,'(2I4,4E14.7)'),i,j, gr(i,j),gv_new(-2:imax,-2:jmax), T_new(-2:imax,-2:jmax)u(i,j),gv(i,j),ge(i,j)!madhu
            !call apply_bc(fr,fu,fv,fe,nx,ny,imax,jmax)
            !call apply_bc(gr,gu,gv,ge,nx,ny,imax,jmax)
        !end do  
      !end do 
        !rho(i,ny) = rho(i,ny-1)*b0+rho(i,1)*bp
        !u(i,ny) =     u(i,ny-1)*b0+  u(i,1)*bp
        !v(i,ny) =     v(i,ny-1)*b0+  v(i,1)*bp
        !T(i,ny) =     T(i,ny-1)*b0+  T(i,1)*bp
        rho(i,ny+1) = rho(i,ny-1)*b0+rho(i,2)*bp
        u(i,ny+1) =     u(i,ny-1)*b0+  u(i,2)*bp
        v(i,ny+1) =     v(i,ny-1)*b0+  v(i,2)*bp
        T(i,ny+1) =     T(i,ny-1)*b0+  T(i,2)*bp
        u(i,ny:ny+3)=1
        v(i,ny:ny+3)=0
        T(i,ny:ny+3)=1
        rho(i,ny:ny+3)=1
        !u(i,ny+0)=0
        !v(i,ny+0)=0
        !u(i,ny+1)=-u(i,ny-1)
        !v(i,ny+1)=-v(i,ny-1)
        !u(i,ny+0)=1
        !u(i,ny+1)=1
        !v(i,ny+0)=0
        !v(i,ny+1)=0
    end do
    !u(0,:)=1.
    !u(nx+1,:)=1.
    !u(:,0)=0
    !u(:,ny+1)=1
    !T(0,:)=T(1,:)
    !T(nx+1,:)=T(nx,:)
    !T(:,0)=T(:,1)
    !T(:,ny+1)=T(:,ny)
    u(:,-2:1)=0
    v(:,-2:1)=0
    T(:,-2:1)=1
    do i=-2,imax
        temp=rho(i,2)
        rho(i,-2)=temp
        rho(i,-1)=temp
        rho(i,0)=temp
        rho(i,1)=temp
    end do
    !u(-2:1,:)=u(nx-3:nx,:)
    !v(-2:1,:)=v(nx-3:nx,:)
    !T(-2:1,:)=T(nx-3:nx,:)
    !rho(-2:1,:)=rho(nx-3:nx,:)
    
    dx=1./(nx-1)
    dy=.1/(ny-1)
    do j = 2, ny-1
        do i = 2, nx-1
          !temp_rho(i,j)=(8.*rho(i,j)+rho(i-1,j)+rho(i+1,j)+rho(i,j-1)+rho(i,j+1))/12.
        end do  
        !u1=u(1:nx,j)
        !v1=v(1:nx,j)
        !call explicitFilter4x(u1, uf, dx, nx, 88, 88)   
        !call explicitFilter4x(v1, vf, dx, nx, 88, 88)
        !u(1:nx,j)=uf 
        !v(1:nx,j)=vf 
    end do
    !rho(-2:imax,-2:jmax)=1.0
    !u=1.0
    !explicitFilter4x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    
end subroutine apply_bc
!----------------------------
subroutine compute_pressure_energy(rho, u, v, T, p, E,re,pr,gamma,mach,s1,nx,ny,imax,jmax)
    implicit none
    double precision, intent(in) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
                                        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
    double precision, intent(out) :: p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
    double precision gamma, mach, s1, re, pr
    integer :: i, j, nx, ny, imax, jmax

    do j = 1, ny
        do i = 1, nx
            p(i,j) = rho(i,j) * T(i,j)/gamma/mach**2
            E(i,j) = p(i,j)/(gamma - 1.0) + 0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2)
            !p(i,j) = (E(i,j)-0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))*(gamma-1)
        end do
    end do
end subroutine compute_pressure_energy
!-------------------------------------
subroutine write_output(x,y,rho, u, v, T, p,nx,ny,imax,jmax,time)
    implicit none
    double precision, intent(in) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
                    v(-2:imax,-2:jmax), T(-2:imax,-2:jmax), p(-2:imax,-2:jmax)
    double precision :: x(-2:nx+3),y(-2:ny+3)
    double precision :: time
    integer :: i, j, nx, ny, imax, jmax
    open(10, file="outfileuv.txt")
    open(11,file="outfile.txt")
    write(10, '(A)') "x,y,u,v,rho,T"
    open(11, file="outfile.txt")
    write(11, '(A)') "x,y,rho,u,v,T"
    do j=1, ny
        do i = 1,nx
           write(10,'(6F15.8)') x(i),y(j),u(i,j),v(i,j),rho(i,j),T(i,j)
        end do
       end do
    do i = 1, nx
        do j = 1, ny
            write(11,'(6E15.7)') x(i), y(j), rho(i,j), u(i,j), v(i,j),  T(i,j)
        end do
    end do
    close(10)
end subroutine write_output
!---------------------------
subroutine inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
    grinv,guinv,gvinv,geinv,nx,ny,imax,jmax)
implicit none
integer :: i, j, nx, ny, imax, jmax
double precision, intent(in) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
double precision, intent(in) ::  p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
double precision :: frinv(-2:imax,-2:jmax), fuinv(-2:imax,-2:jmax),&
grinv(-2:imax,-2:jmax), guinv(-2:imax,-2:jmax)
double precision :: fvinv(-2:imax,-2:jmax), feinv(-2:imax,-2:jmax),&
gvinv(-2:imax,-2:jmax), geinv(-2:imax,-2:jmax)

open(13,file="outfileinvscidflux.txt")

do j = 1, ny
    do i = 1, nx
        frinv(i,j) = (rho(i,j)*u(i,j))
        grinv(i,J) = (rho(i,j)*v(i,j))

        fuinv(i,j) = (rho(i,j)*u(i,j)**2 + p(i,j))
        guinv(i,j) = (rho(i,j)*u(i,j)*v(i,j) + p(i,j))

        fvinv(i,j) = (rho(i,j)*u(i,j)*v(i,j)+ p(i,j))
        gvinv(i,j) = (rho(i,j)*v(i,j)**2 + p(i,j))

        feinv(i,j)=u(i,j)*(E(i,j)+p(i,j))
        geinv(i,j)=v(i,j)*(E(i,j)+p(i,j))
        !write(13,'(2I4,4F16.10)'),i,j, feinv(i,j),fuinv(i,j),fvinv(i,j),feinv(i,j)!madhu
        !write(13,'(2I4,4F16.10)'),i,j, geinv(i,j),guinv(i,j),gvinv(i,j),geinv(i,j)!madhu
    end do
end do

end subroutine inviscid_flux
!---------------------------
subroutine viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
     re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
implicit none
double precision, intent(in) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
double precision, intent(in) ::  p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)

double precision :: gamma, mach, s1, re, pr
integer :: i, j, nx, ny, imax, jmax
double precision, dimension(-2:imax,-2:jmax) :: dudx, dudy, dvdx, dvdy, dTdx, dTdy
double precision :: tauxx, tauxy,tauyy, qx,qy
double precision :: fuvis(-2:imax,-2:jmax),guvis(-2:imax,-2:jmax)
double precision :: fvvis(-2:imax,-2:jmax),gvvis(-2:imax,-2:jmax)
double precision :: fevis(-2:imax,-2:jmax),gevis(-2:imax,-2:jmax)

double precision :: eckpv, amucv
double precision :: div,sxx,syy,sxy
double precision dx,dy,dt
open(14,file="outfileviscousflux.txt")

eckpv=(gamma-1.)*mach*mach*pr

do j = 1, ny
  do i = 1, nx
    !print *,i,j,u(i,j),v(i,j)
    ! Central difference for gradients
    dudx(i,j) = (u(i+1,j) - u(i-1,j)) / (2.0*dx)
    dvdx(i,j) = (v(i+1,j) - v(i-1,j)) / (2.0*dx)
    dTdx(i,j) = (T(i+1,j) - T(i-1,j)) / (2.0*dx)
    dudy(i,j) = (u(i,j+1) - u(i,j-1)) / (2.0*dy)
    dvdy(i,j) = (v(i,j+1) - v(i,j-1)) / (2.0*dy)
    dTdy(i,j) = (T(i,j+1) - T(i,j-1)) / (2.0*dy)
    !dudx(i,j) = (u(i-2,j)-6*u(i-1,j)+3*u(i,j)+2*u(i+1,j))/(6.*dx)
    !dvdx(i,j) = (v(i-2,j)-6*v(i-1,j)+3*v(i,j)+2*v(i+1,j))/(6.*dx)
    !dTdx(i,j) = (T(i-2,j)-6*T(i-1,j)+3*T(i,j)+2*T(i+1,j))/(6.*dx)
    !dudy(i,j) = (u(i,j-2)-6*u(i,j-1)+3*u(i,j)+2*u(i,j+1))/(6.*dy)
    !dvdy(i,j) = (v(i,j-2)-6*v(i,j-1)+3*v(i,j)+2*v(i,j+1))/(6.*dy)
    !dTdy(i,j) = (T(i,j-2)-6*T(i,j-1)+3*T(i,j)+2*T(i,j+1))/(6.*dy)
  end do
  i=1
  dudx(i,j) = (u(i+1,j) - u(i,j)) / (1.0*dx)
  dvdx(i,j) = (v(i+1,j) - v(i,j)) / (1.0*dx)
  dTdx(i,j) = (T(i+1,j) - T(i,j)) / (1.0*dx)
  i=nx
  dudx(i,j) = (u(i,j) - u(i-1,j)) / (1.0*dx)
  dvdx(i,j) = (v(i,j) - v(i-1,j)) / (1.0*dx)
  dTdx(i,j) = (T(i,j) - T(i-1,j)) / (1.0*dx)
end do 
do i = 1, nx
    j=1
    dudy(i,j) = (u(i,j+1) - u(i,j)) / (1.0*dy)
    dvdy(i,j) = (v(i,j+1) - v(i,j)) / (1.0*dy)
    dTdy(i,j) = (T(i,j+1) - T(i,j)) / (1.0*dy)
    j=ny
    dudy(i,j) = (u(i,j) - u(i,j-1)) / (1.0*dy)
    dvdy(i,j) = (v(i,j) - v(i,j-1)) / (1.0*dy)
    dTdy(i,j) = (T(i,j) - T(i,j-1)) / (1.0*dy)
end do 
do j = 1, ny
  do i = 1, nx
    sxx   = (dudx(i,j)+dudx(i,j))
    syy   = (dvdy(i,j)+dvdy(i,j))
    sxy   = (dvdx(i,j)+dudy(i,j))
    div   = (sxx+syy)/2.

! Viscous stress contributions
    amucv=T(i,j)**3./2.*(1.+s1)/(T(i,j)+s1)
    tauxx = amucv/re* (-2./3.*div+sxx)
    tauyy = amucv/re* (-2./3.*div+syy)
    tauxy = amucv/re* (sxy)

! Heat conduction
    qx=-amucv/(gamma-1)/mach**2/re/pr*dTdx(i,j)
    qy=-amucv/(gamma-1)/mach**2/re/pr*dTdy(i,j)

    fuvis(i,j)=tauxx
    guvis(i,j)=tauxy
    fvvis(i,j)=tauxy
    gvvis(i,j)=tauyy
    fevis(i,j)=-qx+u(i,j)*tauxx+v(i,j)*tauxy
    gevis(i,j)=-qy+u(i,j)*tauxy+v(i,j)*tauyy

!write(14,'(2I4,5F13.8)'),i,j, sxy,div,tauxx,tauxy,tauyy !madhu
!write(14,'(2I4,3F16.10)'),i,j, fuvis(i,j),fvvis(i,j),fevis(i,j)!madhu
!write(14,'(2I4,3F16.10)'),i,j, guvis(i,j),gvvis(i,j),gevis(i,j)!madhu
    !if(abs(sxx+syy+sxy+div) > 0.0) then
        !print *,i,j, du_dy,dv_dx,div
        !write(*,'(2I3,6E14.6)') i,j, dudy(i,j),div,u(i,j+1),u(i,j),u(i,j-1),guvis(i,j)
    !continue
    !end if
  end do
end do
end subroutine viscous_flux
!-------------------------------------
subroutine RK_update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new,re,pr,gamma,mach,&
    s1,dt,dx,dy,nx,ny,imax,jmax)
implicit none
integer :: i, j, ir, nx, ny, imax, jmax, nmax54
double precision  rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
double precision   p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
double precision  rho_new(-2:imax,-2:jmax), u_new(-2:imax,-2:jmax), &
v_new(-2:imax,-2:jmax), T_new(-2:imax,-2:jmax)
double precision :: gamma, mach, s1, re, pr
double precision :: ru_new,rv_new,E_new(-2:imax,-2:jmax),p_new,rr_new,rE_new

double precision :: frinv(-2:imax,-2:jmax), frvis(-2:imax,-2:jmax),&
grinv(-2:imax,-2:jmax), grvis(-2:imax,-2:jmax)
double precision :: fuinv(-2:imax,-2:jmax), fuvis(-2:imax,-2:jmax),&
guinv(-2:imax,-2:jmax), guvis(-2:imax,-2:jmax)
double precision :: fvinv(-2:imax,-2:jmax), fvvis(-2:imax,-2:jmax),&
gvinv(-2:imax,-2:jmax), gvvis(-2:imax,-2:jmax)
double precision :: feinv(-2:imax,-2:jmax), fevis(-2:imax,-2:jmax),&
geinv(-2:imax,-2:jmax), gevis(-2:imax,-2:jmax)
double precision :: fu(-2:imax,-2:jmax), fv(-2:imax,-2:jmax),fe(-2:imax,-2:jmax), &
            fr(-2:imax,-2:jmax)
double precision :: gu(-2:imax,-2:jmax), gv(-2:imax,-2:jmax),ge(-2:imax,-2:jmax), &
            gr(-2:imax,-2:jmax)
double precision dx,dy,dt
double precision A54(5),B54(5),C54(5)
    !----------------------------------------------------
    nmax54=5
    A54(1) = 0.0d0
    A54(2) = -567301805773.0d0 / 1357537059087.0d0
    A54(3) = 0.0d0 - 2404267990393.0d0 / 2016746695238.0d0
    A54(4) = 0.0d0 - 3550918686646.0d0 / 2091501179385.0d0
    A54(5) = 0.0d0 - 1275806237668.0d0 / 842570457699.0d0
    B54(1) = 1432997174477.0d0 / 9575080441755.0d0
    B54(2) = 5161836677717.0d0 / 13612068292357.0d0
    B54(3) = 1720146321549.0d0 / 2090206949498.0d0
    B54(4) = 3134564353537.0d0 / 4481467310338.0d0
    B54(5) = 2277821191437.0d0 / 14882151754819.0d0
    C54(2)=1./5.
    C54(3)=3./10.
    C54(4)=3./5.
    C54(1)=0.0
    !----------------------------------------------------

open(12,file="outfiletemp.txt")
rho_new = rho
u_new = u
v_new = v
E_new= E
do ir=1,nmax54
  call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
                          frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
  call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
                  re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
  do j = 1, ny-0
    do i = 1, nx-0
        fr(i,j)=-frinv(i,j)
        fu(i,j)=-fuinv(i,j)+fuvis(i,j)
        fv(i,j)=-fvinv(i,j)+fvvis(i,j)
        fe(i,j)=-feinv(i,j)+fevis(i,j)

        gr(i,j)=-grinv(i,j)
        gu(i,j)=-guinv(i,j)+guvis(i,j)
        gv(i,j)=-gvinv(i,j)+gvvis(i,j)
        ge(i,j)=-geinv(i,j)+gevis(i,j)
    end do  
end do 
 do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        rho_new(i,j)=A54(ir)*rho_new(i,j)          +dt*((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
        u_new(i,j)=A54(ir)*rho_new(i,j)*u_new(i,j) +dt*((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
        v_new(i,j)=A54(ir)*rho_new(i,j)*v_new(i,j) +dt*((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
        E_new(i,j)=A54(ir)*E_new(i,j)              +dt*((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)
 
        rho(i,j)=B54(ir) * rho_new(i,j) + rho(i,j)
        u(i,j)  =B54(ir) * rho_new(i,j)*u_new(i,j)+  rho(i,j)*u(i,j)
        v(i,j)  =B54(ir) * rho_new(i,j)*v_new(i,j)+  rho(i,j)*v(i,j)
        E(i,j)  =B54(ir) * E_new(i,j) + E(i,j)

        p(i,j)=(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T_new(i,j)=p(i,j)/rho_new(i,j)*gamma*mach**2

        u(i,j)=u(i,j)/rho(i,j)
        v(i,j)=v(i,j)/rho(i,j)
        p(i,j)=(E(i,j)-0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2))*(gamma-1)
        T(i,j)=p(i,j)/rho(i,j)*gamma*mach**2

        !write(12,'(2I4,6F16.10)'),i,j,u(i,j),v(i,j),u_new(i,j),v_new(i,j),ru,ru_new
        write(12,'(2I4,5E15.7)'),i,j,rho(i,j),u(i,j),v(i,j),E(i,j),T(i,j)
        !write(12,'(2I4,5E15.7)'),i,j,rho_new(i,j),u_new(i,j),v_new(i,j),E_new(i,j),T_new(i,j)
    end do
  end do
  call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
end do
rho_new = rho
u_new = u
v_new = v
E_new= E
end subroutine RK_update
!--------------------------------------------------------------
subroutine RK2_update(rho, u, v, T, p, E, rho_new, u_new, v_new,T_new,re,pr,gamma,mach,&
    s1,dt,dx,dy,nx,ny,imax,jmax)
implicit none
integer :: i, j, nx, ny, imax, jmax
double precision  rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
double precision  rho0(-2:imax,-2:jmax), u0(-2:imax,-2:jmax), &
        v0(-2:imax,-2:jmax), E0(-2:imax,-2:jmax)
double precision   p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
double precision  rho_new(-2:imax,-2:jmax), u_new(-2:imax,-2:jmax), &
v_new(-2:imax,-2:jmax), T_new(-2:imax,-2:jmax)
double precision :: gamma, mach, s1, re, pr
double precision :: E_new(-2:imax,-2:jmax)

double precision :: frinv(-2:imax,-2:jmax), &
grinv(-2:imax,-2:jmax)
double precision :: fuinv(-2:imax,-2:jmax), fuvis(-2:imax,-2:jmax),&
guinv(-2:imax,-2:jmax), guvis(-2:imax,-2:jmax)
double precision :: fvinv(-2:imax,-2:jmax), fvvis(-2:imax,-2:jmax),&
gvinv(-2:imax,-2:jmax), gvvis(-2:imax,-2:jmax)
double precision :: feinv(-2:imax,-2:jmax), fevis(-2:imax,-2:jmax),&
geinv(-2:imax,-2:jmax), gevis(-2:imax,-2:jmax)
double precision :: fu(-2:imax,-2:jmax), fv(-2:imax,-2:jmax),fe(-2:imax,-2:jmax), &
            fr(-2:imax,-2:jmax)
double precision :: gu(-2:imax,-2:jmax), gv(-2:imax,-2:jmax),ge(-2:imax,-2:jmax), &
            gr(-2:imax,-2:jmax)
double precision :: kr(-2:imax,-2:jmax), ku(-2:imax,-2:jmax),kv(-2:imax,-2:jmax), &
            ke(-2:imax,-2:jmax)
double precision :: mr(-2:imax,-2:jmax), mu(-2:imax,-2:jmax),mv(-2:imax,-2:jmax), &
            me(-2:imax,-2:jmax)            
double precision dx,dy,dt

open(12,file="outfiletemp.txt")
rho_new = rho
u_new = u
v_new = v
E_new= E
rho0 = rho
u0= u
v0 = v
E0= E

  call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
                          frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
  call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
                  re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
  do j = 1, ny-0
    do i = 1, nx-0
        fr(i,j)=-frinv(i,j)
        fu(i,j)=-fuinv(i,j)+fuvis(i,j)
        fv(i,j)=-fvinv(i,j)+fvvis(i,j)
        fe(i,j)=-feinv(i,j)+fevis(i,j)

        gr(i,j)=-grinv(i,j)
        gu(i,j)=-guinv(i,j)+guvis(i,j)
        gv(i,j)=-gvinv(i,j)+gvvis(i,j)
        ge(i,j)=-geinv(i,j)+gevis(i,j)
    end do  
end do 
 do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        kr(i,j)=((fr(i+1,j)-fr(i-0,j))/1./dx+(gr(i,j+1)-gr(i,j-0))/1./dy)
        ku(i,j)=((fu(i+1,j)-fu(i-0,j))/1./dx+(gu(i,j+1)-gu(i,j-0))/1./dy)
        kv(i,j)=((fv(i+1,j)-fv(i-0,j))/1./dx+(gv(i,j+1)-gv(i,j-0))/1./dy)
        ke(i,j)=((fe(i+1,j)-fe(i-0,j))/1./dx+(ge(i,j+1)-ge(i,j-0))/1./dy)
        rho_new(i,j)=rho(i,j)+dt*kr(i,j)
        u_new(i,j)=rho(i,j)*u(i,j)+dt*ku(i,j)
        v_new(i,j)=rho(i,j)*v(i,j)+dt*kv(i,j)
        E_new(i,j)=E(i,j)+dt*ke(i,j)
        rho(i,j)=rho_new(i,j)
        u(i,j)=u_new(i,j)/rho_new(i,j)
        v(i,j)=v_new(i,j)/rho_new(i,j)
        E(i,j)=E_new(i,j)
        p(i,j)=(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T(i,j)=p(i,j)/rho_new(i,j)*gamma*mach**2
    end do
end do
call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny,imax,jmax)
call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
do j = 1, ny-0
 do i = 1, nx-0
  fr(i,j)=-frinv(i,j)
  fu(i,j)=-fuinv(i,j)+fuvis(i,j)
  fv(i,j)=-fvinv(i,j)+fvvis(i,j)
  fe(i,j)=-feinv(i,j)+fevis(i,j)

  gr(i,j)=-grinv(i,j)
  gu(i,j)=-guinv(i,j)+guvis(i,j)
  gv(i,j)=-gvinv(i,j)+gvvis(i,j)
  ge(i,j)=-geinv(i,j)+gevis(i,j)
 end do  
end do 
do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        mr(i,j)=((fr(i+0,j)-fr(i-1,j))/1./dx+(gr(i,j+0)-gr(i,j-1))/1./dy)
        mu(i,j)=((fu(i+0,j)-fu(i-1,j))/1./dx+(gu(i,j+0)-gu(i,j-1))/1./dy)
        mv(i,j)=((fv(i+0,j)-fv(i-1,j))/1./dx+(gv(i,j+0)-gv(i,j-1))/1./dy)
        me(i,j)=((fe(i+0,j)-fe(i-1,j))/1./dx+(ge(i,j+0)-ge(i,j-1))/1./dy)
        rho_new(i,j)=rho0(i,j)      +dt/2.*(kr(i,j)+mr(i,j))
        u_new(i,j)=rho0(i,j)*u0(i,j)+dt/2.*(ku(i,j)+mu(i,j))
        v_new(i,j)=rho0(i,j)*v0(i,j)+dt/2.*(kv(i,j)+mv(i,j))
        E_new(i,j)=E0(i,j)          +dt/2.*(ke(i,j)+me(i,j))
        rho(i,j)=rho_new(i,j)
        u(i,j)=u_new(i,j)/rho_new(i,j)
        v(i,j)=v_new(i,j)/rho_new(i,j)
        E(i,j)=E_new(i,j)
        p(i,j)=(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T(i,j)=p(i,j)/rho_new(i,j)*gamma*mach**2
    end do
end do
call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
!write(12,'(2I4,6F16.10)'),i,j,u(i,j),v(i,j),u_new(i,j),v_new(i,j),ru,ru_new
!write(12,'(2I4,5E15.7)'),i,j,rho(i,j),u(i,j),v(i,j),E(i,j),T(i,j)
!write(12,'(2I4,5E15.7)'),i,j,rho_new(i,j),u_new(i,j),v_new(i,j),E_new(i,j),T_new(i,j)
call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny,imax,jmax)
rho_new = rho
u_new = u
v_new = v
E_new= E
!T_new=T
end subroutine RK2_update
!------------------------------------------
    ! Stage 1
!call rhs_func(u, dudt)
!u1 = u + dt * dudt
! Stage 2
!call rhs_func(u1, dudt)
!u2 = 0.75 * u + 0.25 * (u1 + dt * dudt)
! Stage 3
!call rhs_func(u2, dudt)
!u  = (1.0/3.0) * u + (2.0/3.0) * (u2 + dt * dudt)
!--------------------------------------------------------------
subroutine RK3_update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new,re,pr,gamma,mach,&
    s1,dt,dx,dy,nx,ny,imax,jmax)
implicit none
integer :: i, j, ir, nx, ny, imax, jmax, nmax54
double precision  rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
double precision   p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
double precision  rho_new(-2:imax,-2:jmax), u_new(-2:imax,-2:jmax), &
v_new(-2:imax,-2:jmax), T_new(-2:imax,-2:jmax)
double precision :: gamma, mach, s1, re, pr
double precision :: ru_new,rv_new,E_new(-2:imax,-2:jmax),p_new,rr_new,rE_new

double precision :: frinv(-2:imax,-2:jmax), frvis(-2:imax,-2:jmax),&
grinv(-2:imax,-2:jmax), grvis(-2:imax,-2:jmax)
double precision :: fuinv(-2:imax,-2:jmax), fuvis(-2:imax,-2:jmax),&
guinv(-2:imax,-2:jmax), guvis(-2:imax,-2:jmax)
double precision :: fvinv(-2:imax,-2:jmax), fvvis(-2:imax,-2:jmax),&
gvinv(-2:imax,-2:jmax), gvvis(-2:imax,-2:jmax)
double precision :: feinv(-2:imax,-2:jmax), fevis(-2:imax,-2:jmax),&
geinv(-2:imax,-2:jmax), gevis(-2:imax,-2:jmax)
double precision :: fu(-2:imax,-2:jmax), fv(-2:imax,-2:jmax),fe(-2:imax,-2:jmax), &
            fr(-2:imax,-2:jmax)
double precision :: gu(-2:imax,-2:jmax), gv(-2:imax,-2:jmax),ge(-2:imax,-2:jmax), &
            gr(-2:imax,-2:jmax)
double precision :: kr(-2:imax,-2:jmax), ku(-2:imax,-2:jmax),kv(-2:imax,-2:jmax), &
            ke(-2:imax,-2:jmax)
double precision :: mr(-2:imax,-2:jmax), mu(-2:imax,-2:jmax),mv(-2:imax,-2:jmax), &
            me(-2:imax,-2:jmax)
double precision :: nr(-2:imax,-2:jmax), nu(-2:imax,-2:jmax),nv(-2:imax,-2:jmax), &
            ne(-2:imax,-2:jmax)               
double precision dx,dy,dt
double precision A54(5),B54(5),C54(5)
    !----------------------------------------------------
    nmax54=5
    A54(1) = 0.0d0
    A54(2) = -567301805773.0d0 / 1357537059087.0d0
    A54(3) = 0.0d0 - 2404267990393.0d0 / 2016746695238.0d0
    A54(4) = 0.0d0 - 3550918686646.0d0 / 2091501179385.0d0
    A54(5) = 0.0d0 - 1275806237668.0d0 / 842570457699.0d0
    B54(1) = 1432997174477.0d0 / 9575080441755.0d0
    B54(2) = 5161836677717.0d0 / 13612068292357.0d0
    B54(3) = 1720146321549.0d0 / 2090206949498.0d0
    B54(4) = 3134564353537.0d0 / 4481467310338.0d0
    B54(5) = 2277821191437.0d0 / 14882151754819.0d0
    C54(2)=1./5.
    C54(3)=3./10.
    C54(4)=3./5.
    C54(5)=1.0
    C54(1)=0.0
    !----------------------------------------------------

open(12,file="outfiletemp.txt")
rho_new = rho
u_new = u
v_new = v
E_new= E
!---------------Step 1
  call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
                          frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
  call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
                  re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
  do j = 1, ny-0
    do i = 1, nx-0
        fr(i,j)=-frinv(i,j)
        fu(i,j)=-fuinv(i,j)+fuvis(i,j)
        fv(i,j)=-fvinv(i,j)+fvvis(i,j)
        fe(i,j)=-feinv(i,j)+fevis(i,j)

        gr(i,j)=-grinv(i,j)
        gu(i,j)=-guinv(i,j)+guvis(i,j)
        gv(i,j)=-gvinv(i,j)+gvvis(i,j)
        ge(i,j)=-geinv(i,j)+gevis(i,j)
    end do  
end do 
 do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        kr(i,j)=((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
        ku(i,j)=((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
        kv(i,j)=((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
        ke(i,j)=((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)
        rho_new(i,j)=rho(i,j)       +dt*kr(i,j)
        u_new(i,j)  =rho(i,j)*u(i,j)+dt*ku(i,j)
        v_new(i,j)  =rho(i,j)*v(i,j)+dt*kv(i,j)
        E_new(i,j)  =E(i,j)         +dt*ke(i,j)
        rho(i,j)=rho_new(i,j)
        u(i,j)=u_new(i,j)/rho_new(i,j)
        v(i,j)=v_new(i,j)/rho_new(i,j)
        E(i,j)=E_new(i,j)
        p(i,j)=(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T(i,j)=p(i,j)/rho_new(i,j)*gamma*mach**2
    end do
end do
!------------Step 2
call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny,imax,jmax)
call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
do j = 1, ny-0
do i = 1, nx-0
fr(i,j)=-frinv(i,j)
fu(i,j)=-fuinv(i,j)+fuvis(i,j)
fv(i,j)=-fvinv(i,j)+fvvis(i,j)
fe(i,j)=-feinv(i,j)+fevis(i,j)

gr(i,j)=-grinv(i,j)
gu(i,j)=-guinv(i,j)+guvis(i,j)
gv(i,j)=-gvinv(i,j)+gvvis(i,j)
ge(i,j)=-geinv(i,j)+gevis(i,j)
end do  
end do 
do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        mr(i,j)=((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
        mu(i,j)=((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
        mv(i,j)=((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
        me(i,j)=((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)
        rho_new(i,j)=rho(i,j)       +dt/2.*(kr(i,j)+mr(i,j))/2.
        u_new(i,j)  =rho(i,j)*u(i,j)+dt/2.*(ku(i,j)+mu(i,j))/2.
        v_new(i,j)  =rho(i,j)*v(i,j)+dt/2.*(kv(i,j)+mv(i,j))/2
        E_new(i,j)  =E(i,j)         +dt/2.*(ke(i,j)+me(i,j))/2.
        rho(i,j)=rho_new(i,j)
        u(i,j)  =u_new(i,j)/rho_new(i,j)
        v(i,j)  =v_new(i,j)/rho_new(i,j)
        E(i,j)  =E_new(i,j)
        p(i,j)=(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T(i,j)=p(i,j)/rho_new(i,j)*gamma*mach**2
    end do
end do
!-----------Step 3
call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
call inviscid_flux(rho, u, v, T, p, E,frinv,fuinv,fvinv,feinv,&
frinv,fuinv,fvinv,feinv,nx,ny,imax,jmax)
! Viscous stress contributions
call viscous_flux(rho, u, v, T, p, E, fuvis,fvvis,fevis,guvis,gvvis,gevis, &
re,pr,gamma,mach,s1,dt,dx,dy,nx,ny,imax,jmax)
do j = 1, ny-0
do i = 1, nx-0
fr(i,j)=-frinv(i,j)
fu(i,j)=-fuinv(i,j)+fuvis(i,j)
fv(i,j)=-fvinv(i,j)+fvvis(i,j)
fe(i,j)=-feinv(i,j)+fevis(i,j)

gr(i,j)=-grinv(i,j)
gu(i,j)=-guinv(i,j)+guvis(i,j)
gv(i,j)=-gvinv(i,j)+gvvis(i,j)
ge(i,j)=-geinv(i,j)+gevis(i,j)
end do  
end do 
do j = 2, ny-1
    do i = 2, nx-1
! Update equations
        nr(i,j)=((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
        nu(i,j)=((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
        nv(i,j)=((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
        ne(i,j)=((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)
        rho_new(i,j)=rho(i,j)       +dt/6.*(kr(i,j)+mr(i,j)+4.*nr(i,j))
        u_new(i,j)  =rho(i,j)*u(i,j)+dt/6.*(ku(i,j)+mu(i,j)+4.*nu(i,j))
        v_new(i,j)  =rho(i,j)*v(i,j)+dt/6.*(kv(i,j)+mv(i,j)+4.*nv(i,j))
        E_new(i,j)  =E(i,j)         +dt/6.*(ke(i,j)+me(i,j)+4.*ne(i,j))
        rho(i,j)=rho_new(i,j)
        u(i,j)  =u_new(i,j)/rho_new(i,j)
        v(i,j)  =v_new(i,j)/rho_new(i,j)
        E(i,j)  =E_new(i,j)
        p(i,j)  =(E_new(i,j)-0.5*rho_new(i,j)*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)      
        T(i,j)  =p(i,j)/rho_new(i,j)*gamma*mach**2
    end do
end do
call apply_bc(rho, u, v, T,nx,ny,imax,jmax)
!write(12,'(2I4,6F16.10)'),i,j,u(i,j),v(i,j),u_new(i,j),v_new(i,j),ru,ru_new
!write(12,'(2I4,5E15.7)'),i,j,rho(i,j),u(i,j),v(i,j),E(i,j),T(i,j)
!write(12,'(2I4,5E15.7)'),i,j,rho_new(i,j),u_new(i,j),v_new(i,j),E_new(i,j),T_new(i,j)
call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny,imax,jmax)
rho_new = rho
u_new = u
v_new = v
E_new= E
end subroutine RK3_update
!***************************************************************
subroutine explicitFilter4x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: c0,c1,c2,c3,c4,c5,c6
    double precision :: ft(-5:imax+6), fixed_value_L, fixed_value_R
    integer :: i

    c6 = 0.001254597714
    c5 =-0.008520738659
    c4 = 0.029662754736
    c3 =-0.069975429105
    c2 = 0.123632891797
    c1 =-0.171503832236
    c0 = 0.190899511506

    fixed_value_L = f(1)
    fixed_value_R = f(imax)

    ft(1:imax)=f

    if (boundary_flag_L == 88) then
        ft(0) = f(2)
        ft(-1) = f(3)
        ft(-2) = f(4)
        ft(-3) = f(5)
        ft(-4) = f(6)
        ft(-5) = f(7)       
    end if

    if (boundary_flag_R == 88) then
        ft(imax+1)=f(imax-1)
        ft(imax+2)=f(imax-2)
        ft(imax+3)=f(imax-3)
        ft(imax+4)=f(imax-4)
        ft(imax+5)=f(imax-5)
        ft(imax+6)=f(imax-6)
    end if

    if (boundary_flag_L == 99) then
        ft(-5:0)=f(imax-6:imax-1)
    end if

    if (boundary_flag_R == 99) then
        ft(imax+1:imax+6)=f(2:7)
    end if

    if (boundary_flag_L == 96) then
        ft(1) = 0.0d0
        ft(0) = -f(2)
        ft(-1) = -f(3)
        ft(-2) = -f(4)
        ft(-3) = -f(5)
        ft(-4) = -f(6)
        ft(-5) = -f(7)  
    end if

    if (boundary_flag_R == 96) then
        ft(imax) = 0.0d0
        ft(imax+1)=-f(imax-1)
        ft(imax+2)=-f(imax-2)
        ft(imax+3)=-f(imax-3)
        ft(imax+4)=-f(imax-4)
        ft(imax+5)=-f(imax-5)
        ft(imax+6)=-f(imax-6)
    end if

    if (boundary_flag_L == 100) then
        ft(-5:1) = fixed_value_L
    end if

    if (boundary_flag_R == 100) then
        ft(imax+1:imax+6)=fixed_value_R
    end if

    if (boundary_flag_L < 80) then
        ft(-5:1) = boundary_flag_L
    end if

    if (boundary_flag_R < 80) then
        ft(imax+0:imax+6) = boundary_flag_R
    end if

    do i = 1, imax
        fx(i) = ft(i)+ 0.1*&
               (c6*(ft(i+6)+ft(i-6)) + &
                c5*(ft(i+5)+ft(i-5)) + &
                c4*(ft(i+4)+ft(i-4)) + &
                c3*(ft(i+3)+ft(i-3)) + &
                c2*(ft(i+2)+ft(i-2)) + &
                c1*(ft(i+1)+ft(i-1)) + &
                c0*(ft(i+0)))
    end do

end subroutine explicitFilter4x
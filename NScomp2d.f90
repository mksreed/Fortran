program navier_stokes_fd
    implicit none
    ! Parameters
    integer, parameter :: nx = 20, ny = 20
    integer, parameter :: imax = nx+3, jmax = ny+3
    integer, parameter :: nt = 10000
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax=10, ymax=0.5
    double precision, parameter :: dt = 1e-5
    double precision, parameter :: gamma = 1.4, mach=0.4, s1=110.4/273., re=100, pr=0.72
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

    ! Time-stepping
    time=0
    do n = 1, nt
        call update(rho, u, v, T, p, E,rho_new,u_new,v_new,T_new,re,pr,gamma,mach, &
                                                                s1,dt,dx,dy,nx,ny,imax,jmax)
        call apply_bc(rho_new, u_new, v_new, T_new,nx,ny,imax,jmax)
        call compute_pressure_energy(rho_new,u_new,v_new,T_new,p,E,re,pr,gamma,mach, &
                                                                s1,nx,ny,imax,jmax)

        rho = rho_new
        u = u_new
        v = v_new
        T = T_new
        time=time+dt
        write(*,'(4F12.6)'),time,maxval(u),maxval(v),maxval(p)
    end do

    call write_output(x,y,rho, u, v, T, p,nx,ny,imax,jmax,time)
end program navier_stokes_fd
!---------------------
subroutine initialize(x,y,dx,dy,xmax,ymax,rho, u, v, T,nx,ny,imax,jmax)
    implicit none
    double precision, intent(out) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), v(-2:imax,-2:jmax), &
                                                 T(-2:imax,-2:jmax)
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax, ymax
    double precision :: rad, theta
    integer :: i, j, nx, ny, imax, jmax

    dx=xmax/(nx-1)
    dy=ymax/(ny-1)
    do i=-2,imax,1;x(i)=(i-1)*dx;end do
    do j=-2,jmax,1;y(j)=(j-1)*dy; end do

    do j = -1, ny+1
        do i = -1, nx+1
            rho(i,j) = 1.0
            u(i,j) = 1
            v(i,j) = 0.0
            T(i,j) = 1.0
            theta=atan(y(j),x(i))
            rad=sqrt((x(i))**2+(y(j))**2)
            !u(i,j)=u(i,j)+0.-rad*sin(theta)*exp(-5.*rad)*10
            !v(i,j)=v(i,j)+rad*cos(theta)*exp(-5.*rad)*10
            !write(*,'(2I4,4F20.13)'),i,j,u(i,j),v(i,j),rad,theta
        end do
    end do
    !u=1.
    
end subroutine initialize
!--------------------------------------
subroutine apply_bc(rho, u, v, T,nx,ny,imax,jmax)
    implicit none
    double precision, intent(inout) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
                                        v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
    integer :: i, j, nx, ny, imax, jmax
    integer :: bp=0,b0=1,bs=0,bn=0,ba=0;

    do j = 1, ny
        rho(1,j) = rho(2,j)*b0+rho(nx,j)*bp
        u(1,j)   =   u(2,j)*b0+  u(nx,j)*bp
        v(1,j)   =   v(2,j)*b0+  v(nx,j)*bp
        T(1,j)   =   T(2,j)*b0+  T(nx,j)*bp
        rho(0,j) = rho(2,j)*b0+rho(nx-1,j)*bp
        u(0,j)   =   u(2,j)*b0+  u(nx-1,j)*bp
        v(0,j)   =   v(2,j)*b0+  v(nx-1,j)*bp
        T(0,j)   =   T(2,j)*b0+  T(nx-1,j)*bp
        u(0,j)=1.
        v(0,j)=0.
        u(1,j)=1.
        v(1,j)=0.

        rho(nx,j) = rho(nx-1,j)*b0+rho(1,j)*bp
        u(nx,j)   =   u(nx-1,j)*b0+  u(1,j)*bp
        v(nx,j)   =   v(nx-1,j)*b0+  v(1,j)*bp
        T(nx,j)   =   T(nx-1,j)*b0+  T(1,j)*bp
        rho(nx+1,j) = rho(nx-1,j)*b0+rho(2,j)*bp
        u(nx+1,j)   =   u(nx-1,j)*b0+  u(2,j)*bp
        v(nx+1,j)   =   v(nx-1,j)*b0+  v(2,j)*bp
        T(nx+1,j)   =   T(nx-1,j)*b0+  T(2,j)*bp
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
        u(i,1)=0
        v(i,1)=0
        u(i,0)=-u(i,2)
        v(i,0)=-v(i,2)
        !u(i,0:1)=1

        rho(i,ny) = rho(i,ny-1)*b0+rho(i,1)*bp
        u(i,ny) =     u(i,ny-1)*b0+  u(i,1)*bp
        v(i,ny) =     v(i,ny-1)*b0+  v(i,1)*bp
        T(i,ny) =     T(i,ny-1)*b0+  T(i,1)*bp
        rho(i,ny+1) = rho(i,ny-1)*b0+rho(i,2)*bp
        u(i,ny+1) =     u(i,ny-1)*b0+  u(i,2)*bp
        v(i,ny+1) =     v(i,ny-1)*b0+  v(i,2)*bp
        T(i,ny+1) =     T(i,ny-1)*b0+  T(i,2)*bp
        !u(i,ny+0)=0
        !v(i,ny+0)=0
        !u(i,ny+1)=-u(i,ny-1)
        !v(i,ny+1)=-v(i,ny-1)
        u(i,ny+0)=1
        u(i,ny+1)=1

    end do
    u(0,:)=1.
    u(nx+1,:)=1.
    u(:,0)=0
    u(:,ny+1)=1
    do j = 0, ny+1
        do i = 0, nx+1
          !print *,i,j,u(i,j),v(i,j)
        end do    
    end do
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
            p(i,j) = (E(i,j)-0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))*(gamma-1)
        end do
    end do
end subroutine compute_pressure_energy
!---------------------------
subroutine update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new,re,pr,gamma,mach,&
                                            s1,dt,dx,dy,nx,ny,imax,jmax)
    implicit none
    integer :: i, j, nx, ny, imax, jmax
    double precision, intent(in) :: rho(-2:imax,-2:jmax), u(-2:imax,-2:jmax), &
                                                v(-2:imax,-2:jmax), T(-2:imax,-2:jmax)
    double precision, intent(in) ::  p(-2:imax,-2:jmax), E(-2:imax,-2:jmax)
    double precision, intent(out) :: rho_new(-2:imax,-2:jmax), u_new(-2:imax,-2:jmax), &
                                        v_new(-2:imax,-2:jmax), T_new(-2:imax,-2:jmax)
    double precision :: gamma, mach, s1, re, pr
    double precision :: ru_new,rv_new,e_new,p_new

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
    open(12,file="outfiletemp.txt")

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

            !write(*,'(2I4,4F16.10)'),i,j, fr(i,j),fu(i,j),fv(i,j),fe(i,j)!madhu
            !write(*,'(2I4,4F16.10)'),i,j, gr(i,j),gu(i,j),gv(i,j),ge(i,j)!madhu
            !call apply_bc(fr,fu,fv,fe,nx,ny,imax,jmax)
            !call apply_bc(gr,gu,gv,ge,nx,ny,imax,jmax)

    do j = 2, ny-1
        do i = 2, nx-1
            ! Update equations
            rho_new(i,j) = rho(i,j)+dt*((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
            ru_new=rho(i,j)*u(i,j) +dt*((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
            rv_new=rho(i,j)*v(i,j) +dt*((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
            e_new        =   E(i,j)+dt*((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)

            u_new(i,j) = ru_new/rho_new(i,j)
            v_new(i,j) = rv_new/rho_new(i,j)
            p_new=(e_new-0.5*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)
            T_new(i,j)=p_new/rho_new(i,j)*gamma*mach**2

            !write(12,'(2I4,6F16.10)'),i,j,u(i,j),v(i,j),u_new(i,j),v_new(i,j),ru,ru_new
            !write(12,'(2I4,6F16.10)'),i,j,rho(i,j),u(i,j),v(i,j),E(i,j),T(i,j),p(i,j)
            !write(12,'(2I4,6F16.10)'),i,j,rho_new(i,j),u_new(i,j),v_new(i,j),e_new,T_new(i,j),p_new
        end do
    end do

end subroutine update
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
    write(10, '(A)') "x,y,u,v,t"
    open(11, file="outfile.txt")
    write(11, '(A)') "x,y,r,u,v,p"
    do j=1, ny
        do i = 1,nx
           write(10,'(5F12.6)') x(i),y(j),u(i,j),v(i,j),time
        end do
       end do
    do i = 1, nx
        do j = 1, ny
            write(11,'(6F12.6)') x(i), y(j), rho(i,j), u(i,j), v(i,j),  p(i,j)
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

        fvinv(i,j) = (rho(i,j)*u(i,j+1)*v(i,j)+ p(i,j))
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
double precision :: du_dx, du_dy, dv_dx, dv_dy, dT_dx, dT_dy
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
    du_dx = (u(i+1,j) - u(i-1,j)) / (2.0*dx)
    du_dy = (u(i,j+1) - u(i,j-1)) / (2.0*dy)
    dv_dx = (v(i+1,j) - v(i-1,j)) / (2.0*dx)
    dv_dy = (v(i,j+1) - v(i,j-1)) / (2.0*dy)
    dT_dx = (T(i+1,j) - T(i-1,j)) / (2.0*dx)
    dT_dy = (T(i,j+1) - T(i,j-1)) / (2.0*dy)
    sxx   = (du_dx+du_dx)
    syy   = (dv_dy+dv_dy)
    sxy   = (dv_dx+du_dy)
    div   = (sxx+syy)/2.
if(abs(sxx+syy+sxy+div) > 0.0) then
    !print *,i,j, du_dy,dv_dx,div
    !write(*,'(2I3,5F12.6)') i,j, du_dy,dv_dx,div,u(i,j+1),u(i,j-1)
!continue
end if
! Viscous stress contributions
    amucv=T(i,j)**3./2.*(1.+s1)/(T(i,j)+s1)/re
    tauxx = amucv* (-2./3.*div+sxx)
    tauyy = amucv* (-2./3.*div+syy)
    tauxy = amucv* (sxy)
    fuvis(i,j)=tauxx
    guvis(i,j)=tauxy
    fvvis(i,j)=tauxy
    gvvis(i,j)=tauyy
! Heat conduction
    qx=-amucv/(gamma-1)/mach**2/re/pr*dT_dx
    qy=-amucv/(gamma-1)/mach**2/re/pr*dT_dy
    fevis(i,j)=-qx+u(i,j)*tauxx+v(i,j)*tauxy
    gevis(i,j)=-qy+u(i,j)*tauxy+v(i,j)*tauyy

!write(14,'(2I4,5F13.8)'),i,j, sxy,div,tauxx,tauxy,tauyy !madhu
!write(14,'(2I4,3F16.10)'),i,j, fuvis(i,j),fvvis(i,j),fevis(i,j)!madhu
!write(14,'(2I4,3F16.10)'),i,j, guvis(i,j),gvvis(i,j),gevis(i,j)!madhu

  end do
end do
end subroutine viscous_flux
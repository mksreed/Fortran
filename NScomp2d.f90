program navier_stokes_fd
    implicit none
    ! Parameters
    integer, parameter :: nx = 10, ny = 10
    integer, parameter :: i_max = nx+2, j_max = ny+2
    integer, parameter :: nt = 1000
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax=1, ymax=1
    double precision, parameter :: dt = 1e-5
    double precision, parameter :: gamma = 1.4, mach=0.2, s1=110.4/273., re=100, pr=0.72
    double precision, parameter :: R = 287.0
    double precision, parameter :: mu = 1.8e-5
    double precision, parameter :: kappa = 0.025


    ! State variables
    double precision :: rho(i_max,j_max), u(i_max,j_max), v(i_max,j_max), T(i_max,j_max)
    double precision :: p(i_max,j_max), E(i_max,j_max)
    double precision :: rho_new(i_max,j_max), u_new(i_max,j_max), v_new(i_max,j_max), T_new(i_max,j_max)

    integer :: n

    ! Initialization
    call initialize(x,y,dx,dy,xmax,ymax,rho, u, v, T,nx,ny)
    call compute_pressure_energy(rho, u, v, T, p, E, re, pr, gamma,mach,s1,nx,ny)

    ! Time-stepping
    do n = 1, nt
        call update(rho, u, v, T, p, E,rho_new,u_new,v_new,T_new,re,pr,gamma,mach,s1,dt,dx,dy,nx,ny)
        call apply_bc(rho_new, u_new, v_new, T_new,nx,ny)
        call compute_pressure_energy(rho_new,u_new,v_new,T_new,p,E,re,pr,gamma,mach,s1,nx,ny)

        rho = rho_new
        u = u_new
        v = v_new
        T = T_new
    end do

    call write_output(rho, u, v, T, p,nx,ny)
end program navier_stokes_fd
!---------------------
subroutine initialize(x,y,dx,dy,xmax,ymax,rho, u, v, T,nx,ny)
    implicit none
    double precision, intent(out) :: rho(nx,ny), u(nx,ny), v(nx,ny), T(nx,ny)
    double precision :: x(-2:nx+3),y(-2:ny+3), dx, dy, xmax, ymax
    integer :: i, j, nx, ny

    dx=xmax/(nx-1)
    dy=ymax/(ny-1)
    do i=-2,nx+2,1;x(i)=(i-1)*dx;end do
    do j=-2,ny+2,1;y(j)=(j-1)*dy; end do

    do j = 1, ny
        do i = 1, nx
            rho(i,j) = 1.0
            u(i,j) = 0.0
            v(i,j) = 0.0
            T(i,j) = 1.0
        end do
    end do
end subroutine initialize
!--------------------------------------
subroutine apply_bc(rho, u, v, T,nx,ny)
    implicit none
    double precision, intent(inout) :: rho(nx,ny), u(nx,ny), v(nx,ny), T(nx,ny)
    integer :: i, j, nx, ny

    do j = 1, ny
        rho(1,j) = rho(2,j)
        u(1,j) = u(2,j)
        v(1,j) = v(2,j)
        T(1,j) = T(2,j)

        rho(nx,j) = rho(nx-1,j)
        u(nx,j) = u(nx-1,j)
        v(nx,j) = v(nx-1,j)
        T(nx,j) = T(nx-1,j)
    end do

    do i = 1, nx
        rho(i,1) = rho(i,2)
        u(i,1) = u(i,2)
        v(i,1) = v(i,2)
        T(i,1) = T(i,2)

        rho(i,ny) = rho(i,ny-1)
        u(i,ny) = u(i,ny-1)
        v(i,ny) = v(i,ny-1)
        T(i,ny) = T(i,ny-1)
    end do
end subroutine apply_bc
!----------------------------
subroutine compute_pressure_energy(rho, u, v, T, p, E,re,pr,gamma,mach,s1,nx,ny)
    implicit none
    double precision, intent(in) :: rho(nx,ny), u(nx,ny), v(nx,ny), T(nx,ny)
    double precision, intent(out) :: p(nx,ny), E(nx,ny)
    double precision gamma, mach, s1, re, pr
    integer :: i, j, nx, ny

    do j = 1, ny
        do i = 1, nx
            p(i,j) = rho(i,j) * T(i,j)/gamma/mach**2
            E(i,j) = p(i,j)/(gamma - 1.0) + 0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2)
        end do
    end do
end subroutine compute_pressure_energy
!---------------------------
subroutine update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new,re,pr,gamma,mach,s1,dt,dx,dy,nx,ny)
    implicit none
    double precision, intent(in) :: rho(nx,ny), u(nx,ny), v(nx,ny), T(nx,ny), p(nx,ny), E(nx,ny)
    double precision, intent(out) :: rho_new(nx,ny), u_new(nx,ny), v_new(nx,ny), T_new(nx,ny)
    double precision :: gamma, mach, s1, re, pr
    integer :: i, j, nx, ny
    double precision :: du_dx, du_dy, dv_dx, dv_dy, dT_dx, dT_dy
    double precision :: ru,rv,ruu,rvv,ruv,ru_new,rv_new,e_new,p_new
    double precision :: druu_dx,druv_dx,druv_dy,drvv_dy
    double precision :: tauxx, tauxy,tauyy, qx,qy
    double precision :: frinv(nx,ny), frvis(nx,ny),grinv(nx,ny), grvis(nx,ny)
    double precision :: fuinv(nx,ny), fuvis(nx,ny),guinv(nx,ny), guvis(nx,ny)
    double precision :: fvinv(nx,ny), fvvis(nx,ny),gvinv(nx,ny), gvvis(nx,ny)
    double precision :: feinv(nx,ny), fevis(nx,ny),geinv(nx,ny), gevis(nx,ny)
    double precision :: fu(nx,ny), fv(nx,ny),fe(nx,ny), fr(nx,ny)
    double precision :: gu(nx,ny), gv(nx,ny),ge(nx,ny), gr(nx,ny)
    double precision :: eckpv, amucv

    double precision :: s11,s22,s12,div,sxx,syy,sxy

    double precision dx,dy,dt

    eckpv=(gamma-1.)*mach*mach*pr
    do j = 2, ny-1
        do i = 2, nx-1
            ru = rho(i,j)*u(i,j)
            rv = rho(i,j)*v(i,j)
            ruu = rho(i,j)*u(i,j)*u(i,j)
            rvv = rho(i,j)*v(i,j)*v(i,j)
            ruv = rho(i,j)*u(i,j)*v(i,j)
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
            div   = (s11+s22)/2.
            if(abs(sxx+syy+sxy+div) > 0.0) then
                print *,i,j, sxy,div
                !continue
            end if
            frinv(i,j) = (rho(i,j)*u(i,j))
            grinv(i,J) = (rho(i,j)*v(i,j))

            fuinv(i,j) = (rho(i,j)*u(i,j)**2 + p(i,j))
            guinv(i,j) = (rho(i,j)*u(i,j)*v(i,j) + p(i,j))

            fvinv(i,j) = (rho(i,j)*u(i,j+1)*v(i,j)+ p(i,j))
            gvinv(i,j) = (rho(i,j)*v(i,j)**2 + p(i,j))


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

            feinv(i,j)=u(i,j)*E(i,j)+p(i,j)
            geinv(i,j)=v(i,j)*E(i,j)+p(i,j)
            fevis(i,j)=-qx+u(i,j)*tauxx+v(i,j)*tauxy
            gevis(i,j)=-qy+u(i,j)*tauxy+v(i,j)*tauyy


            fr(i,j)=-frinv(i,j)
            fu(i,j)=-fuinv(i,j)+fuvis(i,j)
            fv(i,j)=-fvinv(i,j)+fvvis(i,j)
            fe(i,j)=-feinv(i,j)+fevis(i,j)

            gr(i,j)=-grinv(i,j)
            gu(i,j)=-guinv(i,j)+guvis(i,j)
            gv(i,j)=-gvinv(i,j)+gvvis(i,j)
            ge(i,j)=-geinv(i,j)+gevis(i,j)
            !write(*,'(2I4,5F16.10)'),i,j, sxy,div,tauxx,tauxy,tauyy !madhu
            !write(*,'(2I4,4F16.10)'),i,j, fuinv(i,j),fuvis(i,j),guinv(i,j),guvis(i,j)!madhu
            !write(*,'(2I4,4F16.10)'),i,j, fvinv(i,j),fvvis(i,j),gvinv(i,j),gvvis(i,j)!madhu
            !write(*,'(2I4,4F16.10)'),i,j, feinv(i,j),fevis(i,j),geinv(i,j),gevis(i,j)!madhu
            !write(*,'(2I4,4F16.10)'),i,j, fr(i,j),fu(i,j),fv(i,j),fe(i,j)!madhu
            !write(*,'(2I4,4F16.10)'),i,j, gr(i,j),gu(i,j),gv(i,j),ge(i,j)!madhu
            call apply_bc(fr,fu,fv,fe,nx,ny)
            call apply_bc(gr,gu,gv,ge,nx,ny)
        end do
    end do
    do j = 2, ny-1
        do i = 2, nx-1
            ! Update equations
            rho_new(i,j) = rho(i,j)+dt*((fr(i+1,j)-fr(i-1,j))/2./dx+(gr(i,j+1)-gr(i,j-1))/2./dy)
            ru_new       =       ru+dt*((fu(i+1,j)-fu(i-1,j))/2./dx+(gu(i,j+1)-gu(i,j-1))/2./dy)
            rv_new       =       rv+dt*((fv(i+1,j)-fv(i-1,j))/2./dx+(gv(i,j+1)-gv(i,j-1))/2./dy)
            e_new        =   E(i,j)+dt*((fe(i+1,j)-fe(i-1,j))/2./dx+(ge(i,j+1)-ge(i,j-1))/2./dy)
            !p(i,j) = rho(i,j) * T(i,j)/gamma/mach**2
            u_new(i,j) = ru_new/rho_new(i,j)
            v_new(i,j) = rv_new/rho_new(i,j)
            p_new=(e_new-0.5*(u_new(i,j)**2+v_new(i,j)**2))*(gamma-1)
            T_new(i,j)=p_new/rho_new(i,j)*gamma*mach**2
            !write(*,'(2I4,6F16.10)'),i,j,rho(i,j),u(i,j),v(i,j),E(i,j),T(i,j),p(i,j)
            !write(*,'(2I4,6F16.10)'),i,j,rho_new(i,j),u_new(i,j),v_new(i,j),e_new,T_new(i,j),p_new
        end do
    end do
end subroutine update
!-------------------------------------
subroutine write_output(rho, u, v, T, p,nx,ny)
    implicit none
    double precision, intent(in) :: rho(nx,ny), u(nx,ny), v(nx,ny), T(nx,ny), p(nx,ny)
    integer :: i, j, nx, ny
    open(10, file="solution.dat")
    do j = 2, ny-1
        do i = 2, nx-1
            write(10,*) i, j, rho(i,j), u(i,j), v(i,j), T(i,j), p(i,j)
        end do
    end do
    close(10)
end subroutine write_output
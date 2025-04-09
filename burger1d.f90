program call_thomas_algorithm
    implicit none
    integer :: i,ir, n, nmax,nmax54,imax,i_print_freq
    double precision :: pi, xmax, dx, s1, s2, s3,nu,dt,time
    double precision :: A54(5), B54(5), C54(5)! RK 5 
    integer, parameter :: nx = 65  ! Define nx as needed
    double precision :: x(-1:nx+2)
    double precision :: uf(nx)
    double precision :: u1(nx),u0(nx),fdudx(nx)
    double precision :: u(-1:nx+2),ue(-1:nx+2)
    double precision :: alpha=1./3, aarhs=14./9.,bbrhs=1./9.
    integer :: boundary_flag_L, boundary_flag_R
    integer :: unit=10, unit1=20
    character(len=20) :: filename, filename1
  !1 0                          1432997174477/9575080441755
  !2 567301805773/1357537059087 5161836677717/13612068292357
 !3 2404267990393/2016746695238 1720146321549/2090206949498
  !4 3550918686646/2091501179385 3134564353537/4481467310338
 !5 1275806237668/842570457699 2277821191437/14882151754819
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
    pi = 3.14159265358979323d0
    !-----------------------------------------------------
    nu=0.001
    nmax=8001
    time=0.0
    dt=0.0001
    i_print_freq=nmax/5
    imax=nx
    xmax = 2.0d0 * pi
    dx = xmax / real(nx - 1)
    s1 = 2.0d0
    s2=1.0
    s3=1.0
    !x=[(i*dx-dx,i=-1,nx+2,1)]
    x=[((i-nx/2)*dx,i=-1,nx+2,1)]
    u=max(cos(s1*x)+cos(s2*x)*s3,0.0)
    ue=u

    u0=u(1:nx)
    u1=u0
    filename = "outfile.txt"
    filename1 = "outfile1.txt"
    open(unit=unit, file=filename, status='replace')
    open(unit=unit1, file=filename1, status='replace')

    boundary_flag_L=99
    boundary_flag_R=99

    !--------------------------------------------------------
        ! Main loop
    !u0=0
    !u1=0
    uf=0
    write(unit1, '(A)') "t,umax"
    write(unit, '(A)') "i,x,u,ue,t"
    do n = 1, nmax + 1
        do ir = 1, nmax54
            fdudx=0
            call compute_R(u0, time, fdudx, dx, imax, nu, boundary_flag_L, boundary_flag_R)
            !fdudx=time+C54(ir)*dt
            do i = 1, imax 
                u1(i) = A54(ir) * u1(i) + dt * fdudx(i)
            end do
            do i = 0, imax
                u0(i) = u0(i) + B54(ir) * u1(i)
            end do
        end do
        !uf=uf+dt*time
        uf=u0
        call explicitFilterx(u0, uf, dx, nx, boundary_flag_L, boundary_flag_R)
        time = time + dt
        u0=uf
        !ue=time*time/2.
        if (mod(n, i_print_freq) == 0 .or. n==1) then
            write(unit1, '(2F14.7)') real(n * dt + 0.000),maxval(u0)
            do i = 1,imax
                write(unit,'(I3,4F12.6)') i, x(i),u0(i),ue(i),time
            end do
        end if
    end do
    u(1:imax)=u0
    !--------------------------------------------------------

    do i = 1, nx
        write(*,'(I3,4F12.6)')i, x(i),u(i),ue(i),time
        write(unit,'(I3,4F12.6)') i, x(i),u(i),ue(i),time
    end do
    !call explicitFilterx(u0, uf, dx, nx, boundary_flag_L, boundary_flag_R)

    !write(unit1, '(A)') "i,x,u,ue"
    do i = 1, nx
        !write(*,'(I3,3F12.6)') i, x(i),u0(i),uf(i)
        !write(unit,'(I3,3F12.6)') i, x(i),u0(i),uf(i)
    end do
    close(unit)
    close(unit1)
end program call_thomas_algorithm
  !*************************************************
subroutine compute_R(ut, time, fdudx, dx, imax, nu, boundary_flag_L, boundary_flag_R)
    implicit none
    double precision, dimension(imax), intent(in) :: ut
    double precision, dimension(imax), intent(out) :: fdudx
    double precision, dimension(imax) :: dudx, dudxx
    double precision :: dx, time
    integer :: i,imax
    double precision :: nu
    integer :: boundary_flag_L, boundary_flag_R

    ! Call the explicit functions
    dudx=0
    call explicit6x(ut, dudx, dx, imax, boundary_flag_L, boundary_flag_R)
    dudxx=0
    call explicit6xx(ut, dudxx, dx, imax, boundary_flag_L, boundary_flag_R)

    ! Compute fdudx 
    do i = 1, imax
        fdudx(i) = -ut(i) * dudx(i) + dudxx(i) * nu
        ! fdudx(i) = C * dudx(i) + dudxx6(i) * nu  ! Uncomment if needed
        !fdudx(i)=time
    end do
end subroutine compute_R
  !*************************************************
  subroutine explicit6x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: cm3, cm2, cm1, c0
    double precision :: cp3, cp2, cp1
    double precision :: ft(-2:imax+3), fixed_value_L, fixed_value_R
    integer :: i

    cm3 = -1.0d0 / 60.0d0
    cm2 = 3.0d0 / 20.0d0
    cm1 = -3.0d0 / 4.0d0
    c0 = 0.0d0
    cp3 = -cm3
    cp2 = -cm2
    cp1 = -cm1

    fixed_value_L = f(1)
    fixed_value_R = f(imax)

    ft(1:imax)=f

    if (boundary_flag_L == 88) then
        ft(0) = f(2)
        ft(-1) = f(3)
        ft(-2) = f(4)
    end if

    if (boundary_flag_R == 88) then
        ft(imax+1)=f(imax-1)
        ft(imax+2)=f(imax-2)
        ft(imax+3)=f(imax-3)
    end if

    if (boundary_flag_L == 99) then
        ft(-2:0)=f(imax-3:imax-1)
    end if

    if (boundary_flag_R == 99) then
        ft(imax+1:imax+3)=f(2:4)
    end if

    if (boundary_flag_L == 96) then
        ft(1) = 0.0d0
        ft(0) = -f(2)
        ft(-1) = -f(3)
        ft(-2) = -f(4)
    end if

    if (boundary_flag_R == 96) then
        ft(imax) = 0.0d0
        ft(imax+1)=-f(imax-1)
        ft(imax+2)=-f(imax-2)
        ft(imax+3)=-f(imax-3)
    end if

    if (boundary_flag_L == 100) then
        ft(-2:1) = fixed_value_L
    end if

    if (boundary_flag_R == 100) then
        ft(imax+1:imax+3)=fixed_value_R
    end if

    if (boundary_flag_L < 80) then
        ft(-2:1) = boundary_flag_L
    end if

    if (boundary_flag_R < 80) then
        ft(imax+0:imax+3) = boundary_flag_R
    end if

    do i = 1, imax
        fx(i) = (cm3 * ft(i - 3) + cm2 * ft(i - 2) + cm1 * ft(i - 1) + c0 * ft(i) + &
                     cp1 * ft(i + 1) + cp2 * ft(i + 2) + cp3 * ft(i + 3)) / dx
    end do

end subroutine explicit6x
!***************************************************************
subroutine explicitFilterx(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: ft(-2:imax+3), fixed_value_L, fixed_value_R
    double precision :: a40, a41, a42, a43, alpha
    integer :: i

    alpha = 0.0
    a40 = 11.0 / 16.0 + 5.0 * alpha / 8.0
    a41 = 15.0 / 32.0 + 17.0 * alpha / 16.0
    a42 = -3.0 / 16.0 + 3.0 * alpha / 8.0
    a43 = 1.0 / 32.0 - 1.0 * alpha / 16.0

    fixed_value_L = f(1)
    fixed_value_R = f(imax)

        ft(1:imax)=f

    if (boundary_flag_L == 88) then
        ft(0) = f(2)
        ft(-1) = f(3)
        ft(-2) = f(4)
    end if

    if (boundary_flag_R == 88) then
        ft(imax+1)=f(imax-1)
        ft(imax+2)=f(imax-2)
        ft(imax+3)=f(imax-3)
    end if

    if (boundary_flag_L == 99) then
        ft(-2:0)=f(imax-3:imax-1)
    end if

    if (boundary_flag_R == 99) then
        ft(imax+1:imax+3)=f(2:4)
    end if

    if (boundary_flag_L == 96) then
        ft(1) = 0.0d0
        ft(0) = -f(2)
        ft(-1) = -f(3)
        ft(-2) = -f(4)
    end if

    if (boundary_flag_R == 96) then
        ft(imax) = 0.0d0
        ft(imax+1)=-f(imax-1)
        ft(imax+2)=-f(imax-2)
        ft(imax+3)=-f(imax-3)
    end if

    if (boundary_flag_L == 100) then
        ft(-2:1) = fixed_value_L
    end if

    if (boundary_flag_R == 100) then
        ft(imax+1:imax+3)=fixed_value_R
    end if

    if (boundary_flag_L < 80) then
        ft(-2:1) = boundary_flag_L
    end if

    if (boundary_flag_R < 80) then
        ft(imax+0:imax+3) = boundary_flag_R
    end if

    do i = 1, imax
         fx(i) =    (a40 * (ft(i) + ft(i)) + &
                     a41 * (ft(i + 1) + ft(i - 1)) + &
                     a42 * (ft(i + 2) + ft(i - 2)) + &
                     a43 * (ft(i + 3) + ft(i - 3))) * 0.5
    end do

end subroutine explicitFilterx
  !*************************************************
subroutine explicit6xx(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: cm3, cm2, cm1, c0
    double precision :: cp3, cp2, cp1
    double precision :: ft(-2:imax+3), fixed_value_L, fixed_value_R
    integer :: i
    !double cm3 = 1. / 90., cm2 = -3. / 20., cm1 = 3. / 2., c0 = -49. / 18;

    cm3 = 1.0d0 / 90.0d0
    cm2 = -3.0d0 / 20.0d0
    cm1 = 3.0d0 / 2.0d0
    c0 = -49.0d0 /18;
    cp3 = cm3
    cp2 = cm2
    cp1 = cm1

    fixed_value_L = f(1)
    fixed_value_R = f(imax)

    ft(1:imax)=f

    if (boundary_flag_L == 88) then
        ft(0) = f(2)
        ft(-1) = f(3)
        ft(-2) = f(4)
    end if

    if (boundary_flag_R == 88) then
        ft(imax+1)=f(imax-1)
        ft(imax+2)=f(imax-2)
        ft(imax+3)=f(imax-3)
    end if

    if (boundary_flag_L == 99) then
        ft(-2:0)=f(imax-3:imax-1)
    end if

    if (boundary_flag_R == 99) then
        ft(imax+1:imax+3)=f(2:4)
    end if

    if (boundary_flag_L == 96) then
        ft(1) = 0.0d0
        ft(0) = -f(2)
        ft(-1) = -f(3)
        ft(-2) = -f(4)
    end if

    if (boundary_flag_R == 96) then
        ft(imax) = 0.0d0
        ft(imax+1)=-f(imax-1)
        ft(imax+2)=-f(imax-2)
        ft(imax+3)=-f(imax-3)
    end if

    if (boundary_flag_L == 100) then
        ft(-2:1) = fixed_value_L
    end if

    if (boundary_flag_R == 100) then
        ft(imax+1:imax+3)=fixed_value_R
    end if

    if (boundary_flag_L < 80) then
        ft(-2:1) = boundary_flag_L
    end if

    if (boundary_flag_R < 80) then
        ft(imax+0:imax+3) = boundary_flag_R
    end if

    do i = 1, imax
        fx(i) = (cm3 * ft(i - 3) + cm2 * ft(i - 2) + cm1 * ft(i - 1) + c0 * ft(i) + &
            cp1 * ft(i + 1) + cp2 * ft(i + 2) + cp3 * ft(i + 3)) / dx/dx
    end do
end subroutine explicit6xx
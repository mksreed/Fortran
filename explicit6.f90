program call_thomas_algorithm
    implicit none
    integer :: i
    double precision :: pi, xmax, dx, s1, s2, s3
    integer, parameter :: nx = 100  ! Define nx as needed
    double precision :: x(-1:nx+2)
    double precision :: uxf(nx),uf(nx),uxxf(nx)
    double precision :: uxxe(-1:nx+2),uxx(nx),uu(nx)
    double precision :: u(-1:nx+2),uxe(-1:nx+2),ux(nx)
    double precision :: alpha=1./3, aarhs=14./9.,bbrhs=1./9.
    integer :: boundary_flag_L, boundary_flag_R
    integer :: unit=10, unit1=20
    character(len=20) :: filename, filename1

    pi = 3.14159265358979323d0
    xmax = 2.0d0 * pi
    dx = xmax / real(nx - 1)
    s1 = 5.0d0
    s2=10.0
    s3=1.0
    x=[(i*dx-dx,i=-1,nx+2,1)]
    u=sin(s1*x)+sin(s2*x)*s3
    uxe=s1*cos(s1*x)+s2*cos(s2*x)*s3
    uxxe=-s1*s1*sin(s1*x)-s2*s2*sin(s2*x)*s3
    uu=u(1:nx)
    filename = "outfile.txt"
    filename1 = "outfile1.txt"
    open(unit=unit, file=filename, status='replace')
    open(unit=unit1, file=filename1, status='replace')

    boundary_flag_L=96
    boundary_flag_R=96
    ux=0
    call explicit6x(uu, ux, dx, nx, boundary_flag_L, boundary_flag_R)
    uxx=0
    call explicit6xx(uu, uxx, dx, nx, boundary_flag_L, boundary_flag_R)

    write(unit, '(A)') "i,x,u,ux,uxe,uxx,uxxe"
    do i = 1, nx
        write(*,'(I3,6F12.6)')i, x(i),u(i),ux(i),uxe(i),uxx(i),uxxe(i)
        write(unit,'(I3,6F12.6)') i, x(i),u(i),ux(i),uxe(i),uxx(i),uxxe(i)
    end do
    call explicitFilterx(uxx, uxxf, dx, nx, boundary_flag_L, boundary_flag_R)
    write(unit1, '(A)') "i,x,u,uf"
    do i = 1, nx
        write(*,'(I3,3F12.6)') i, x(i),uxx(i),uxxf(i)
        write(unit1,'(I3,3F12.6)') i,x(i),uxx(i),uxxf(i)
    end do
    close(unit)
    close(unit1)
end program call_thomas_algorithm

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
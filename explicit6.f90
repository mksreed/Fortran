program call_thomas_algorithm
    implicit none
    integer :: i
    double precision :: pi, xmax, dx, s1, s2, s3
    integer, parameter :: nx = 65  ! Define nx as needed
    double precision :: x(-1:nx+2)
    double precision :: u4f(nx),u6f(nx),u10f(nx)
    double precision :: uxf(nx),uf(nx),uxxf(nx)
    double precision :: uxxe(-1:nx+2),uxx(nx),uu(nx)
    double precision :: u(-1:nx+2),uxe(-1:nx+2),ux(nx),u4x(nx)
    double precision :: alpha=1./3, aarhs=14./9.,bbrhs=1./9.
    integer :: boundary_flag_L, boundary_flag_R
    integer :: unit=10, unit1=20
    character(len=20) :: filename, filename1

    pi = 3.14159265358979323d0
    xmax = 2.0d0 * pi
    dx = xmax / real(nx - 1)
    s1 = 1.0d0
    s2=8.0
    s3=0.5
    x=[(i*dx-dx,i=-1,nx+2,1)]
    u=sin(s1*x)+sin(s2*x)*s3
    uxe= s1*cos(s1*x)+s2*cos(s2*x)*s3
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
    u4x=0
    call explicit4x(uu, u4x, dx, nx, boundary_flag_L, boundary_flag_R)   
    uxx=0
    call explicit6xx(uu, uxx, dx, nx, boundary_flag_L, boundary_flag_R)
    write(unit, '(A)') "i,x,u,u6x,u4x,uxe,u6xx,uxxe"
    do i = 1, nx
        write(*,'(I3,7F12.6)')i, x(i),u(i),ux(i),u4x(i),uxe(i),uxx(i),uxxe(i)
        write(unit,'(I3,7F12.6)') i, x(i),u(i),ux(i),u4x(i),uxe(i),uxx(i),uxxe(i)
    end do
    call explicitFilter6x(uu, u6f, dx, nx, boundary_flag_L, boundary_flag_R)
    call explicitFilter10x(uu, u10f, dx, nx, boundary_flag_L, boundary_flag_R)
    call explicitFilter4x(uu, u4f, dx, nx, boundary_flag_L, boundary_flag_R)
    boundary_flag_L=100
    write(unit1, '(A)') "i,x,uu,u4f,u6f,u10f"
    do i = 1, nx
        write(*,'(I3,5F12.6)') i, x(i),uu(i),u4f(i),u6f(i),u10f(i)
        write(unit1,'(I3,5F12.6)') i, x(i),uu(i),u4f(i),u6f(i),u10f(i)
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
subroutine explicitFilter6x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
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

end subroutine explicitFilter6x
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
!*************************************************
subroutine explicit4x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: c0,c1,c2,c3,c4,c5,c6
    double precision :: ft(-5:imax+6), fixed_value_L, fixed_value_R
    integer :: i

    c6=-0.000957455525961
    c5= 0.008242459236975
    c4=-0.037162191039544
    c3 =0.119465303396051
    c2 =-0.320910877852970
    c1 = 0.896607046646854
    c0 = 0.0d0

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
        fx(i) = (c6*(ft(i+6)-ft(i-6)) + &
                c5*(ft(i+5)-ft(i-5)) + &
                c4*(ft(i+4)-ft(i-4)) + &
                c3*(ft(i+3)-ft(i-3)) + &
                c2*(ft(i+2)-ft(i-2)) + &
                c1*(ft(i+1)-ft(i-1)) + &
                c0*(ft(i+0)))/dx
    end do

end subroutine explicit4x
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
!***************************************************************
subroutine explicitFilter10x(f, fx, dx, imax, boundary_flag_L, boundary_flag_R)
    implicit none
    integer, intent(in) :: imax, boundary_flag_L, boundary_flag_R
    double precision, intent(in) :: dx
    double precision, dimension(imax), intent(in) :: f
    double precision, dimension(imax), intent(out) :: fx
    double precision :: c0,c1,c2,c3,c4,c5,c6
    double precision :: ft(-5:imax+6), fixed_value_L, fixed_value_R
    integer :: i

    c6 = 0.0
    c5 = 1./512.
    c4 =-5./256.
    c3 = 45./512.
    c2 =-15./64.
    c1 = 105./256.
    c0 = 193./256.

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
        fx(i) =(c6*(ft(i+6)+ft(i-6)) + &
                c5*(ft(i+5)+ft(i-5)) + &
                c4*(ft(i+4)+ft(i-4)) + &
                c3*(ft(i+3)+ft(i-3)) + &
                c2*(ft(i+2)+ft(i-2)) + &
                c1*(ft(i+1)+ft(i-1)) + &
                c0*(ft(i+0)+ft(i+0)))/2.
    end do

end subroutine explicitFilter10x
!***************************************************************
program call_thomas_algorithm
    implicit none
    integer :: i
    double precision :: pi, xmax, dx, s1
    integer, parameter :: nx = 100  ! Define nx as needed
    double precision :: a(nx), b(nx), c(nx), x(-1:nx+2), d(nx)
    double precision :: ap(nx), bp(nx), cp(nx),dp(nx),uexf(nx)
    double precision :: aas(nx), bas(nx), cas(nx), das(nx)
    double precision :: ss1(nx), ss2(nx), u(-1:nx+2),ue(-1:nx+2),uexp(nx)
    double precision :: alpha=1./3, aarhs=14./9.,bbrhs=1./9.
    integer :: boundary_flag_L, boundary_flag_R
    integer :: unit=10, unit1=20
    character(len=20) :: filename, filename1

    pi = 3.14159265358979323d0
    xmax = 2.0d0 * pi
    dx = xmax / real(nx - 1)
    s1 = 5.0d0
    x=[(i*dx-dx,i=-1,nx+2,1)]
    u=sin(s1*x)
    ue=s1*cos(s1*x)
    filename = "outfile.txt"
    filename1 = "outfile1.txt"
    open(unit=unit, file=filename, status='replace')
    open(unit=unit1, file=filename1, status='replace')

    print *, "Started thomas algorithm"
    do i = 1, nx
        a(i) = alpha
        b(i) = 1.0
        c(i) = alpha
        d(i) = aarhs/2./dx*(u(i+1)-u(i-1))+bbrhs/4./dx*(u(i+2)-u(i-2))
        ap(i) = a(i)
        bp(i) = b(i)
        cp(i) = c(i)
        dp(i) = d(i)
        aas(i) = a(i)
        bas(i) = b(i)
        cas(i) = c(i)
        das(i) = d(i)
    end do

    i = 1
    a(i) = 0.0d0
    b(i) = 1.0d0
    c(i) = 0.0d0
    d(i) = s1*1.0
    cas(i) = cas(i) + aas(i) !96
    !cas(i) = cas(i) - aas(i) !88

    i = nx
    aas(i) = aas(i) + cas(i) !96
    !aas(i) = aas(i) - cas(i) !88
    a(i) = 0.0d0
    b(i) = 1.0d0
    c(i) = 0.0d0
    d(i) = s1*1.0

    call cyclic_thomas(nx - 0, dp, ap, bp, cp, ss1, ss2)
    !call solve_periodic_tridiagonal(ap, bp, cp, dp, nx, ss1)
    dp(nx) = dp(1)
    !dp=ss1

    call solve_tridiag(a,b,c,d,ss1,nx)
    d=ss1
    ss1=0;
    call solve_tridiag(aas,bas,cas,das,ss1,nx)
    das=ss1

    boundary_flag_L=96
    boundary_flag_R=96

    uexf=u(1:nx)
    call explicit6x(uexf, uexp, dx, nx, boundary_flag_L, boundary_flag_R)

    write(unit, '(A)') "i,x,u,up,ue,uas,uexp"
    do i = 1, nx
        write(*,'(I3,6F12.6)')i, x(i),d(i),dp(i),ue(i),das(i),uexp(i)
        write(unit,'(I3,6F12.6)') i,x(i),d(i),dp(i),ue(i),das(i),uexp(i)
    end do
    call explicitFilterx(uexf, uexp, dx, nx, boundary_flag_L, boundary_flag_R)
    write(unit1, '(A)') "i,x,u,uf"
    do i = 1, nx
        write(*,'(I3,3F12.6)') i, x(i),uexf(i),uexp(i)
        write(unit1,'(I3,3F12.6)') i,x(i),uexf(i),uexp(i)
    end do
    close(unit)
end program call_thomas_algorithm
!*************************************
subroutine cyclic_thomas(nx, x, a, b, c, tempcmod, v)
        integer, intent(in) :: nx
        double precision, intent(inout) :: x(nx)
        double precision, intent(in) :: a(nx), b(nx), c(nx)
        double precision, intent(out) :: tempcmod(nx), v(nx)
        double precision :: m
        integer :: ix

        ! first solve a system of length nx - 1 for two right hand sides, ignoring ix == 0
        tempcmod(2) = c(2) / b(2)
        v(2) = -a(2) / b(2)
        x(2) = x(2) / b(2)

        ! loop from 2 to nx - 1 inclusive
        do ix = 3, nx - 1
            m = 1.0d0 / (b(ix) - a(ix) * tempcmod(ix - 1))
            tempcmod(ix) = c(ix) * m
            v(ix) = (0.0d0 - a(ix) * v(ix - 1)) * m
            x(ix) = (x(ix) - a(ix) * x(ix - 1)) * m
        end do

        ! handle nx - 1
        m = 1.0d0 / (b(nx) - a(nx) * tempcmod(nx - 1))
        tempcmod(nx) = c(nx) * m
        v(nx) = (-c(1) - a(nx) * v(nx - 1)) * m
        x(nx) = (x(nx) - a(nx) * x(nx - 1)) * m

        ! loop from nx - 2 to 1 inclusive
        do ix = nx - 1, 2, -1
            v(ix) = v(ix) - tempcmod(ix) * v(ix + 1)
            x(ix) = x(ix) - tempcmod(ix) * x(ix + 1)
        end do

        x(1) = (x(1) - a(1) * x(nx) - c(1) * x(2)) / (b(1) + a(1) * v(nx) + c(1) * v(2))

        ! loop from 1 to nx - 0 inclusive
        do ix = 2, nx - 0
            x(ix) = x(ix) + x(1) * v(ix)
        end do
end subroutine cyclic_thomas
!--------------------------------------------------------------------------------------
! sub routine for triadiagonalsolver
!--------------------------------------------------------------------------------------
subroutine solve_tridiag(a,b,c,d,x,n)
    implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations

      integer,parameter :: r8 = kind(1.d0)

      integer,intent(in) :: n
      double precision,dimension(n),intent(in) :: a,b,c,d
      double precision,dimension(n),intent(out) :: x
      double precision,dimension(n) :: cp,dp
      double precision :: m
      integer i

! initialize c-prime and d-prime
      cp(1) = c(1)/b(1)
      dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
       do i = 2,n
         m = b(i)-cp(i-1)*a(i)
         cp(i) = c(i)/m
         dp(i) = (d(i)-dp(i-1)*a(i))/m
       end do
! initialize x
       x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
      do i = n-1, 1, -1
        x(i) = dp(i)-cp(i)*x(i+1)
      end do

end subroutine solve_tridiag
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
    double precision :: cm3, cm2, cm1, c0
    double precision :: cp3, cp2, cp1
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
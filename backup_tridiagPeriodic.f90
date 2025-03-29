program call_thomas_algorithm
    implicit none
    integer :: i
    double precision :: pi, xmax, dx, s1
    integer, parameter :: nx = 100  ! Define nx as needed
    double precision :: a(nx), b(nx), c(nx), x(nx), d(nx)
    double precision :: ap(nx), bp(nx), cp(nx), dp(nx)
    double precision :: aas(nx), bas(nx), cas(nx), das(nx)
    double precision :: ss1(nx), ss2(nx)
    integer :: unit
    character(len=20) :: filename

    pi = 3.14159265358979323d0
    xmax = 2.0d0 * pi
    dx = xmax / real(nx - 1)
    s1 = 3.0d0
    filename = "outfile.txt"
    open(unit=unit, file=filename, status='replace')

    print *, "Started thomas algorithm"
    do i = 1, nx
        x(i) = real(i-1) * dx
        a(i) = 1.0d0 / (dx * dx)
        b(i) = -2.0d0 / (dx * dx)
        c(i) = 1.0d0 / (dx * dx)
        d(i) = -s1 * s1 * sin(s1 * x(i))
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
    !d(i) = 1.0d0
    cas(i) = cas(i) - aas(i)


    i = nx
    aas(i) = aas(i) - cas(i)
    a(i) = 0.0d0
    b(i) = 1.0d0
    c(i) = 0.0d0
    !d(i) = 1.0d0

    call cyclic_thomas(nx - 0, dp, ap, bp, cp, ss1, ss2)
    !dp(nx) = dp(1)
    !call solve_periodic_tridiagonal(ap, bp, cp, dp, nx, ss1)
    !dp=ss1
    


    !call thomas(nx, d, a, b, c, ss1)
    call solve_tridiag(a,b,c,d,ss1,nx)
    d=ss1
    call thomas(nx, das, aas, bas, cas, ss1)

    write(unit, '(A)') "i,x,u,up,ue,uas"
    do i = 1, nx
        write(*, '(I3,5F12.6)')i, x(i), d(i), dp(i), sin(s1 * x(i)), das(i)
        write(unit, '(I3,5F12.6)') i, x(i), d(i), dp(i), sin(s1 * x(i)), das(i)
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
! *******************************************************************************
subroutine thomas(nx, x, a, b, c, scratch)
    implicit none
    integer ix
    integer, intent(in) :: nx
    double precision, intent(inout) :: x(nx)
    double precision, intent(in) :: a(nx), b(nx), c(nx)
    double precision, intent(out) :: scratch(nx)
    
    ! solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
    ! nx = number of equations
    ! x() = initially contains the input v, and returns x. indexed from [1, ..., nx]
    ! a() = subdiagonal, indexed from [2, ..., nx]
    ! b() = main diagonal, indexed from [1, ..., nx]
    ! c() = superdiagonal, indexed from [1, ..., nx-1]
    ! scratch() = scratch space of length nx, provided by caller, allowing a, b, c to be const

    scratch(1) = c(1) / b(1)
    x(1) = x(1) / b(1)

    ! loop from 2 to nx inclusive
    do ix = 2, nx
        if (ix < nx) then
            scratch(ix) = c(ix) / (b(ix) - a(ix - 1) * scratch(ix - 1))
        end if
        x(ix) = (x(ix) - a(ix - 1) * x(ix - 1)) / (b(ix) - a(ix - 1) * scratch(ix - 1))
    end do

    ! loop from nx - 1 to 1 inclusive
    do ix = nx - 1, 1, -1
        x(ix) = x(ix) - scratch(ix) * x(ix + 1)
    end do
end subroutine thomas
! *******************************************************************************
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
!*****************************************************************
subroutine solve_periodic_tridiagonal(a, b, c, d, n, x)
    double precision, dimension(n), intent(in) :: a, b, c, d
    integer, intent(in) :: n
    double precision, dimension(n), intent(out) :: x
    double precision, dimension(n) :: alpha, beta
    integer :: i

    ! Initialize alpha and beta
    alpha(1) = c(1) / b(1)
    beta(1) = d(1) / b(1)
    
    ! Forward sweep for Thomas algorithm with periodic boundary
    do i = 2, n
      alpha(i) = c(i) / (b(i) - a(i) * alpha(i-1))
      beta(i) = (d(i) - a(i) * beta(i-1)) / (b(i) - a(i) * alpha(i-1))
    end do
    
    ! Periodic correction for last row
    alpha(n) = c(n) / (b(n) - a(n) * alpha(n-1))
    beta(n) = (d(n) - a(n) * beta(n-1)) / (b(n) - a(n) * alpha(n-1))
    
    ! Back substitution with periodic boundary
    x(n) = (d(n) - a(n) * beta(n)) / (b(n) - a(n) * alpha(n-1))
    do i = n-1, 1, -1
      x(i) = beta(i) - alpha(i) * x(i+1)
    end do
    
  end subroutine solve_periodic_tridiagonal
!
!end program periodic_tridiagonal_solver

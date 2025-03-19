program poisson_main
    integer nx,ny,i;
    nx=10
    ny=10
    do i=1,5
        call poisson(nx,ny)
        nx=nx*2
        ny=ny*2
    end do
end program poisson_main
Subroutine poisson(nx, ny)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, parameter :: nt = 200000
    real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
    real(8), parameter :: ymin = 0.0d0, ymax = 3.141592653589793d0
    real(8) :: dx, dy
    real(8), dimension(nx, ny) :: p, pd, b
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(8), dimension(nx, ny) :: Xs, Ys
    integer :: it, cx, cy,i,j
    real(8), dimension(nx,ny) :: pexact, err, totalerr
    real(8) :: meanerr, maxerr, maxp, maxpexact
    ! Initialization
    dx = (xmax - xmin) / real(nx - 1, 8)
    dy = (ymax - ymin) / real(ny - 1, 8)
    p = 0.0d0
    pd = 0.0d0
    b = 0.0d0
    x = [(xmin + i * dx, i = 0, nx-1)]
    y = [(ymin + j * dy, j = 0, ny-1)]
    do 10 j=1,ny,1
        do 10 i=1,nx,1
            Xs(i,j)=x(i);
            Ys(i,j)=y(j);
010 continue
    cx = 8
    cy = 8
    b(ny / 4 + 1, nx / 4 + 1) = 100.0d0
    b(3 * ny / 4 + 1, 3 * nx / 4 + 1) = -100.0d0
    b = 2.0d0 * ((Ys * Ys - Ys) + (Xs * Xs - Xs))
    pexact = (Xs - Xs * Xs) * (Ys - Ys * Ys)
    pexact = sin(cx * Xs) * sin(cy * Ys)
    b = -1.0d0 * (cx * cx * pexact + cy * cy * pexact)
    do it = 1, nt
        pd = p
        p(2:nx-1, 2:ny-1) = (((pd(2:nx-1, 3:ny) + pd(2:nx-1, 1:ny-2)) * dx**2 + &
                               (pd(3:nx, 2:ny-1) + pd(1:nx-2, 2:ny-1)) * dy**2 - &
                               b(2:nx-1, 2:ny-1) * dx**2 * dy**2) / &
                               (2.0d0 * (dx**2 + dy**2)))
        p(1, :) = 0.0d0
        p(nx, :) = 0.0d0
        p(:, 1) = 0.0d0
        p(:, ny) = 0.0d0
        err = p - pexact
        totalerr = sum(err)
        meanerr = sum(err) / real(ny * nx, 8)
        maxerr = maxval(err)
        maxp = maxval(p)
        maxpexact = maxval(pexact)
    end do
    print *, maxerr, meanerr, maxp, maxpexact
end SUBROUTINE poisson
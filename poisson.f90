program wave_equation
    implicit none
    integer, parameter :: nx = 40, ny = 40, nt = 200000
    real(8), parameter :: xmin = 0.0d0, xmax = 3.141592653589793d0
    real(8), parameter :: ymin = 0.0d0, ymax = 3.141592653589793d0
    real(8) :: dx, dy
    real(8), dimension(nx, ny) :: p1, pd1, b1
    real(8), dimension(0:nx+1,0:ny+1) :: p, pd, b
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(8), dimension(0:nx+1, 0:ny+1) :: Xs, Ys
    integer :: it, cx, cy,i,j
    real(8), dimension(0:nx+1,0:ny+1) :: pexact, err, totalerr
    real(8) :: meanerr, maxerr, maxp, maxpexact

    ! Initialization
    dx = (xmax - xmin) / real(nx - 1, 8)
    dy = (ymax - ymin) / real(ny - 1, 8)

    p = 0.0d0
    pd = 0.0d0
    b = 0.0d0
    x = [(xmin + i * dx, i = 0, nx-1)]
    y = [(ymin + j * dy, j = 0, ny-1)]
    do  j=1,ny,1
        do i=1,nx,1
            Xs(i,j)=x(i);
            Ys(i,j)=y(j);
        end do
    end do


    ! Source
    cx = 8
    cy = 8
    b(ny / 4 + 1, nx / 4 + 1) = 100.0d0
    b(3 * ny / 4 + 1, 3 * nx / 4 + 1) = -100.0d0
    b = 2.0d0 * ((Ys * Ys - Ys) + (Xs * Xs - Xs))
    pexact = (Xs - Xs * Xs) * (Ys - Ys * Ys)
    pexact = sin(cx * Xs) * sin(cy * Ys)
    b = -1.0d0 * (cx * cx * pexact + cy * cy * pexact)
!f_xx = (-1*f[i-2]+16*f[i-1]-30*f[i+0]+16*f[i+1]-1*f[i+2])/(12*1
    
    
    h.0*h**2)
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
end program wave_equation


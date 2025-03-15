    program hello 
        implicit none
        integer i,j
        integer, parameter :: nx = 20, ny = 20, nz=12, neq=20
        real(8) x(0:nx+1), y(ny),z(0:nx+1,ny,nz),f(0:nx+1),fx(nx),fxx(nx)
        real(8) , parameter :: dx=2*4.*atan(1.0)/(nx-1)
        x=[((i-1)*dx,i=0,nx+1,1)]
        f=sin(x)
        fxx=(f(2:nx+1)-2*f(1:nx)+f(0:nx-1))/dx/dx
        fx=(f(2:nx+1)-f(0:nx-1))/2/dx
        z(:,1,3)=x
        y=0;
        y(2:ny-1)=x(1:nx-2)
        print '(f7.3)',fx-cos(x(1:nx)),fxx+sin(x(1:nx))
    end program hello 

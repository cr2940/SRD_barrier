subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Parameters for problem
    real(kind=8), parameter :: a = 1.d0
    real(kind=8), parameter :: sigma = 0.5d0
    real(kind=8), parameter :: h0 = 0.1d0

    ! Other storage
    integer :: i,j,ilo,jlo
    real(kind=8) :: omega,x,y,eta,sl,jump,xp,xm,yp,ym

    omega = sqrt(2.d0 * grav * h0) / a

    ilo = floor((.05d0*dx)/dx)
    jlo = floor((.05d0*dy)/dy)

    sl = 1.2d0
    jump = 0.8d0 !.1.5d0
    q(1,:,:) = sl
    q(2,:,:) = 0.d0
    q(3,:,:) = 0.d0

    do i=1-mbc,mx+mbc
        ! x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
            ! y = ylower + (j - 0.5d0) * dy
            ! eta = sigma * h0 / a**2 * (2.d0 * x - sigma)

            ym = ylower + (jlo+j-1.d0) * dy
            yp = ylower + (jlo+j) * dy
            y = 0.5d0*(ym+yp)
            xm = xlower + (ilo+i-1.d0) * dx
            xp = xlower + (ilo+i) * dx
            x = 0.5d0*(xm+xp)

            ! q(1,i,j) = 1.2d0 !max(0.d0,eta - aux(1,i,j))
            ! q(2,i,j) = 0.d0
            ! q(3,i,j) = 0.d0!sigma * omega * q(1,i,j)

            if (y<0.1) then
              q(1,i,j) = q(1,i,j) + jump
            end if
            if (0.4d0*(x) + 0.299d0 < y .and. &
                   0.4d0*x + 0.301d0 > y ) then
                   q(1,i,j) = 0.d0
                 end if

        enddo
    enddo




end subroutine qinit

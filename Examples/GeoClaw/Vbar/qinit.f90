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
    real(kind=8) :: omega,x,y,eta,ym,yp,xm,xp,sl,jump

    omega = sqrt(2.d0 * grav * h0) / a

    ilo = floor((.05d0*dx)/dx)
    jlo = floor((.05d0*dy)/dy)

    ! sea surface level, at rest:
    sl = 1.2d0
    jump = 0.8d0
    q(1,:,:) = sl
    q(2,:,:) = 0.d0
    q(3,:,:) = 0.d0

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc

          ym = ylower + (jlo+j-1.d0) * dy
          yp = ylower + (jlo+j) * dy
          y = 0.5d0*(ym+yp)
          xm = xlower + (ilo+i-1.d0) * dx
          xp = xlower + (ilo+i) * dx
          x = 0.5d0*(xm+xp)

    ! the dam break jump
          if (y>0.9) then
            q(1,i,j) = q(1,i,j) + jump
          end if
    ! dry condition on top of dam
          if (xp <0.5) then
              if  (-0.614*x + 0.718 < y .and.&   ! 0.72
                 -0.614*x + 0.728> y) then
                      q(1,i,j) = 0.d0
                          end if
            end if
            if (xp > 0.5) then
                        if (0.614*x + 0.104 < y .and. &  ! 0.106
                         0.614*x + 0.114> y) then
                             q(1,i,j) = 0.d0
                          end if
                        end if
        enddo
    enddo


end subroutine qinit

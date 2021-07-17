! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
  ! HLLE solver for the 2D shallow water equations.

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x-momentum
!       3 y-momentum

! Auxiliary quantities:
!       1 b bathymetry
!       2 a_x
!       3 a_y   where (a_x,a_y) is unit normal to left face
!       4 length_ratio_left ratio of length of left face to dyc
!       5 b_x
!       6 b_y   where (b_x,b_y) is unit normal to bottom face
!       7 length_ratio_bottom   ratio of length of bottom face to dxc
!       8 cell_area   ratio of cell area to dxc*dyc
!         (approximately Jacobian of mapping function)
!
!
! Rotate velocity and then call standard Riemann solver.
! The resulting waves and flux differences are then rotated
! back to x-y.

! solve Riemann problems along one slice of data.
!
! On input, ql contains the state vector at the left edge of each cell
! qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
! or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
! f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!     amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.
!
! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                   and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr.
    use redistribute2D ! upload redistribute2d module in the riemann folder
    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq

    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat, v_l, v_r, v_hat, h_hat
    double precision :: h_l, h_r, hsqrt_l, hsqrt_r, hsq2
    double precision :: grav
    double precision :: h_m, hu_m, hv_m, s1, s2
    integer :: depth, mu, mv
    integer :: i, m, mw
    double precision :: unorl, unorr, utanl,utanr, a1, a2,a3
    double precision :: alpha, beta, bar_y,wall_height,wave_wall(3,3),s_wall(3)
    double precision :: amdq_wall(3), apdq_wall(3),delta(3)
    integer ::  inx, iny, ilenrat, mlo
    logical :: L2R, R2L

    common /cparam/  grav

!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv to the orthogonal
!     # momentum:
!
 ! barrier location :
    bar_y = 0.72d0
    wall_height = 1.5d0

    depth = 1
    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    if (ixy.eq.1) then
        inx = 2
        iny = 3
        ilenrat = 4
    else
        inx = 4
        iny = 5
        ilenrat = 7
    endif

! if in column direction, need to do redistribute wave at y=0.72:
     if (ixy == 2) then
       mlo = bar_y * mx
     end if


! Determine rotation matrix:
!               [ alpha  beta ]
!               [-beta  alpha ]
! Note that this reduces to the identity on a standard cartesian grid.
  ! initialize
    wave(:,:,:) = 0.d0
! Determine normal velocity components at this edge:
    do i=2-mbc,mx+mbc
        ! Height
        h_l = qr(depth,i-1)
        h_r = ql(depth,i)

        ! rotation vectors:
        alpha = auxl(inx,i)
        beta  = auxl(iny,i)

        ! rotated Velocity
        unorl = alpha*qr(mu,i-1) + beta*qr(mv,i-1)
        unorr = alpha*ql(mu,i) + beta*ql(mv,i)
        utanl = -beta*qr(mu,i-1) + alpha*qr(mv,i-1)
        utanr = -beta*ql(mu,i) + alpha*ql(mv,i)

        ! delta(1) = h_r-h_l
        ! delta(2) = unorr-unorl
        ! delta(3) = utanr-utanl

        ! EXCUSE me, but if you are now at the barrier edge, i.e. ixy=2 && i=mlo, then you need to
        ! do redistribute_fwave, which will give you wave(:,:,i). Feed the rotated values in, when you
        ! get back out the wave(:,:,i), rotate that as seen below.
        if (ixy == 2 .and. i == mlo) then
          call redistribute_fwave(ixy,(/h_l,utanl,unorl/),(/h_r,utanr,unorr/),&
             auxr(1,i-1),auxl(1,i),wall_height,&
             1,wave_wall,s_wall,amdq_wall,apdq_wall,3,&
             3,L2R,R2L)
        !   ! s(1,i) = s_wall(1) *auxl(ilenrat,i)
        !   ! s(2,i) = s_wall(2) * auxl(ilenrat,i)
        !   ! s(3,i) =s_wall(3) * auxl(ilenrat,i)
        !   ! wave(mu,1,i) = alpha*wave_wall(mu,1) - beta*wave_wall(mv,1)
        !   ! wave(mv,1,i) = beta*wave_wall(mu,1) + alpha*wave_wall(mv,1)
        !   ! ! if (s_wall(2) < 0 ) then
        !   ! ! else
        !   !   ! s(1,i) = s_wall(1)
        !   !   ! s(2,i) = 0.5d0*(s_wall(2)+s_wall(3))
        !   !   wave(mu,2,i) = alpha*wave_wall(mu,2) - beta*wave_wall(mv,2)
        !   !   wave(mv,2,i) = beta*wave_wall(mu,2) + alpha*wave_wall(mv,2)
        !   ! ! end if
        !   ! wave(mu,3,i) = (alpha*wave_wall(mu,3) - beta*wave_wall(mv,3))
        !   ! wave(mv,3,i) = (beta*wave_wall(mu,3) + alpha*wave_wall(mv,3))
          amdq(1,i) = amdq_wall(1)
          amdq(mu,i) = alpha*amdq_wall(mu) - beta*amdq_wall(mv)
          amdq(mv,i) = beta*amdq_wall(mu) + alpha*amdq_wall(mv)
          apdq(1,i) = apdq_wall(1)
          apdq(mu,i) = alpha*apdq_wall(mu) - beta*apdq_wall(mv)
          apdq(mv,i) = beta*apdq_wall(mu) + alpha*apdq_wall(mv)
        !
        else



        ! gravity wave speed
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)

        hsqrt_l = dsqrt(qr(depth,i-1))
        hsqrt_r = dsqrt(ql(depth,i))
        hsq2 = hsqrt_l + hsqrt_r
        h_hat = 0.5*(h_l + h_r)
        u_hat = (unorl/h_l*hsqrt_l + unorr/h_r*hsqrt_r) / hsq2
        v_hat = (utanl/h_l*hsqrt_l + utanr/h_r*hsqrt_r) / hsq2
        c_hat = dsqrt(grav*h_hat)

        ! a3 = -(u_hat-c_hat) * delta(1) + delta(2)
        ! a3 = a3* (0.5d0/c_hat)
        ! a2 = -v_hat * delta(1) + delta(3)
        ! a1 = (u_hat+c_hat)*delta(1) - delta(2)
        ! a1 = a1*(0.5d0/c_hat)
        !
        ! wave(depth,1,i) = a1
        ! wave(mu,1,i) = alpha*a1(u_hat-c_hat) - beta*a1*v_hat
        ! wave(mv,1,i) = beta*a1(u_hat-c_hat) + alpha*a1*v_hat
        ! s(1,i) = (u_hat-c_hat) *  auxl(ilenrat,i)
        !
        ! wave(depth,2,i) = 0.0d0
        ! wave(mu,2,i) = -beta*a2
        ! wave(mv,2,i) = alpha*a2
        ! s(2,i) = u_hat * auxl(ilenrat,i)
        !
        ! wave(1,3,i) = a3
        ! wave(mu,3,i) = alpha*a3*(u_hat+c_hat) - beta*a3*v_hat
        ! wave(mv,3,i) = beta*a3*(u_hat+c_hat) + alpha*a3*v_hat
        ! s(3,i) = (u_hat+c_hat) * auxl(ilenrat,i)


        ! ! Speeds of non-shear waves
        s1 = min(unorl/h_l - c_l, u_hat - c_hat) * auxl(ilenrat,i)
        s2 = max(unorr/h_r + c_r, u_hat + c_hat) *  auxl(ilenrat,i)

        ! "middle" state
        h_m = (unorr - unorl - s2*ql(depth,i) + s1*qr(depth,i-1))/(s1-s2)
        hu_m = (unorr*(unorr/h_r-s2) - unorl*(unorl/h_l-s1) + 0.5*grav*(h_r**2 - h_l**2) ) / (s1-s2)
        hv_m = (utanr*unorr/h_r - utanl*unorl/h_l - s2*utanr + s1*utanl)/(s1-s2)

        wave(depth,1,i) = h_m - h_l
        wave(mu,1,i) = alpha*(hu_m - unorl) - beta*( hv_m - utanl)
        wave(mv,1,i) = beta* (hu_m - unorl) + alpha*(hv_m - utanl)
        s(1,i) = s1

        s(2,i) = u_hat * auxl(ilenrat,i)
        wave(depth,2,i) =0.0d0
        wave(mu,2,i) = -beta*(v_hat*(h_r-h_l) + utanr-utanl)
        wave(mv,2,i) = alpha*(v_hat*(h_r-h_l) + utanr-utanl)


        wave(depth,3,i) = h_r - h_m
        wave(mu,3,i) = alpha*(unorr - hu_m) - beta*(utanr - hv_m)
        wave(mv,3,i) = beta*(unorr - hu_m) + alpha*(utanr - hv_m)
        s(3,i) = s2

      endif

    end do

111 continue
    do m=1, meqn
        do i=2-mbc, mx+mbc
          if (i==mlo) then
            cycle
          end if
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
            end do
        end do
    end do
      ! if (ixy==2) then
      !     print*, "WAVE:", amdq(1,:)
      !   end if

end subroutine rpn2

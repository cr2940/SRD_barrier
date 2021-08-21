module redistribute2D

      implicit none
      real(kind=8):: g,drytol
      ! real(8),parameter :: grav=1.d0
      real(8),parameter :: dry_tolerance=1.e-8
contains

        subroutine barrier_passing(hL,hR,huL,huR,bL,bR,wall_height,&
                  L2R,R2L,hstarL,hstarR,ustarL,ustarR)
            ! determine whether water will overtop a barrier from L2R or R2L
            ! and what intermed height and velocity will be upon impact

            implicit none
            logical L2R, R2L, rare1, rare2
            real(kind=8) :: hL,hR,huL,huR,bL,bR,wall_height,hm,s1m,s2m
            real(kind=8) :: uL,um,uR
            real(kind=8) :: hstarL,hstarR,ustarL,ustarR,hstartest
            real(kind=8) :: ustartest,grav 

            common /cparam/  grav
            g = grav
            drytol = dry_tolerance
            ! print *, "hL,hR: " , hL,hR,g,drytol
            L2R = .false.
            R2L = .false.
            hstarL = 0.d0
            hstarR = 0.d0
            ustarL = 0.d0
            ustarR = 0.d0
            if ( hL.gt.drytol ) then
              uL = huL/hL
              call riemanntype(hL,hL,uL,-uL,hm,um,s1m,s2m,rare1,rare2,&
                                                     5,drytol,g)
              hstartest = max(hL,hm)
              ustartest = max(uL,um)
              ! print *, "HM,BL,BR:",hm,bL,bR
              ! print*, "rise1:",hstartest+bL,min(bL,bR)+wall_height
              if ( hstartest+bL.gt.min(bL,bR)+wall_height ) then
                L2R = .true.
                hstarL = hstartest + bL - wall_height - min(bL,bR)
              end if
              ustarL = ustartest
            end if

            if ( hR.gt.drytol ) then
              uR = huR/hR
              call riemanntype(hR,hR,-uR,uR,hm,um,s1m,s2m,rare1,rare2,&
                                                    5,drytol,g)
              hstartest = max(hR,hm)
              ustartest = max(uR,um)
              ! print*, "rise2:",hstartest+bR,min(bL,bR)+wall_height

              if ( hstartest+bR.gt.min(bL,bR)+wall_height ) then
                R2L = .true.
                hstarR = hstartest + bR - wall_height - min(bL,bR)
              end if
              ustarR = ustartest
            end if
            ! print *, "RESULT OVERTOP:", L2R, R2L
            return
            end subroutine barrier_passing

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine redistribute_fwave(ixy,ql,qr,auxl,auxr,wall_height,&
             maxiter,wave_wall,s_wall,amdq_wall,apdq_wall,mwaves,&
             meqn,L2R,R2L) !turn on for general angled problems
          ! redistributes waves at the barrier to left going and right going fluctuations using augmented solver
          implicit none
          real(kind=8) :: fwave(3,3,2),s(3,2),amdq(3,2),apdq(3,2)
          real(kind=8) :: q_wall(3,3),aux_wall(3),s_wall(3)
          real(kind=8) :: amdq_wall(3),apdq_wall(3),wave_wall(3,3),pL
          real(kind=8) :: ql(3),qr(3),auxl,auxr,wall_height,hvL,hvR,pR
          real(kind=8) :: hstarL,hstarR,ustarL,ustarR,diff1,diff2,hL,hR
          real(kind=8) :: huL,huR,bL,bR,delb,hstar,um,hstartest,phiL,phiR
          real(kind=8) :: uL,uR,vL,vR,wall(3),sL,sR,uhat,chat,sRoe1,sRoe2
          real(kind=8) :: sE1,sE2,swr(mwaves),fwr(meqn,mwaves),s1m,s2m
          real(kind=8) :: del(3),r(3,3),A(3,3),det1,det2,det3,v_av,v_L,v_R
          real(kind=8) :: determinant,beta(3),hustar_interface
          real(kind=8) :: q_wall_l(3,2),q_wall_r(3,2),aux_wall_l(2)
          real(kind=8) :: aux_wall_r(2),grav
          integer :: maxiter,ixy,mu,nv,i,mw,meqn,mwaves,k
          logical :: L2R,R2L,rare1,rare2
          common /cparam/  grav
          ! print *, "ql: ",ql
          ! print *, "qr: ", qr
          ! print *, "WALL HT:", wall_height
          ! print *, "L2R,R2L:", L2R,R2L
          g = grav
          drytol = dry_tolerance

          ! set up orientation
          if (ixy==1) then
             mu=2 ! "x-axis" is normal direction
             nv=3 ! transverse
          else if (ixy ==2 ) then
             mu=3 ! "y-axis" is normal direction
             nv=2
          end if

          ! avoid negative values
          if (qr(1).lt.0.d0) then
                qr(1)=0.d0
                qr(2)=0.d0
                qr(3)=0.d0
          endif

          if (ql(1).lt.0.d0) then
                ql(1)=0.d0
                ql(2)=0.d0
                ql(3)=0.d0
          endif

          ! set up wall problem
          q_wall(1:3,1) = ql
          q_wall(:,2) = 0.d0
          q_wall(1:3,3) = qr

          aux_wall(1) = auxl
          aux_wall(3) = auxr
          aux_wall(2) = min(auxl,auxr) + wall_height

          if ( ql(1) .gt.drytol ) then
            v_L = ql(nv)/ql(1)
          else
            v_L = 0.d0
          end if

          if (qr(1) .gt. drytol ) then
            v_R = qr(nv)/qr(1)
          else
            v_R = 0.d0
          end if

          if (ql(1).lt.drytol .and. qr(1).lt.drytol) then
            v_av = 0.5d0 * (v_L+v_R)
          else
          v_av=(sqrt(g*ql(1))*v_L + sqrt(g*qr(1))*v_R) &
          /(sqrt(g*qr(1))+sqrt(g*ql(1)))  !ROE AV
          end if

          ! overtop or not
          ! ghost state setup
          if ( L2R .and. R2L  ) then
            ! if wall is shorter than the bathymetry on the right
            if ( aux_wall(2) .le. aux_wall(3) .or. wall_height.eq.0.d0) then
              q_wall(1:3,2) = q_wall(1:3,3)
              aux_wall(2) = aux_wall(3)

            ! if wall is shorter " " " " left
          else if ( aux_wall(2) .le. aux_wall(1).or.wall_height.eq.0.d0 ) then
              q_wall(1:3,2) = q_wall(1:3,1)
              aux_wall(2) = aux_wall(1)

            ! wall is higher than both bathymetries
         elseif (aux_wall(2).gt.aux_wall(1).and. &
              aux_wall(2).gt.aux_wall(3)) then
              diff2=aux_wall(2)-aux_wall(3)
              diff1=aux_wall(2)-aux_wall(1)
              q_wall(1,2)=min(q_wall(1,3)-diff2,q_wall(1,1)-diff1)
              q_wall(mu,2)=min(q_wall(mu,3),q_wall(mu,1))
              q_wall(nv,2)=min(q_wall(nv,3),q_wall(nv,1))!0.5d0*(q_wall(nv,3)+q_wall(nv,1))
          end if
        end if
        ! print *,"HTS:", q_wall(1,:)

          ! initialize amdq and apdq
          amdq(1:3,:) = 0.d0
          apdq(1:3,:) = 0.d0
          q_wall_l = q_wall(:,1:2)
          q_wall_r = q_wall(:,2:3)
          aux_wall_l = aux_wall(1:2)
          aux_wall_r = aux_wall(2:3)
          ! riemann problems
          do i = 1, 2
            hL = q_wall_l(1,i)
            hR = q_wall_r(1,i)
            huL = q_wall_l(mu,i)
            huR = q_wall_r(mu,i)
            bL = aux_wall_l(i)
            bR = aux_wall_r(i)
            delb = bR-bL
            hvL = q_wall_l(nv,i)
            hvR = q_wall_r(nv,i)

            !check for wet/dry boundary
            if (hR.gt.drytol) then
               uR=huR/hR
               vR=hvR/hR
               phiR = 0.5d0*g*hR**2 + huR**2/hR
            else
               hR = 0.d0
               huR = 0.d0
               hvR = 0.d0
               uR = 0.d0
               vR = 0.d0
               phiR = 0.d0
            endif

            if (hL.gt.drytol) then
               uL=huL/hL
               vL=hvL/hL
               phiL = 0.5d0*g*hL**2 + huL**2/hL
            else
               hL=0.d0
               huL=0.d0
               hvL=0.d0
               uL=0.d0
               vL=0.d0
               phiL = 0.d0
            endif
            ! throw away the problem if there are no water at all
            if (hL <= drytol .and. hR <= drytol) then
              fwave(:,:,i) = 0.d0
              s(:,i) = 0.d0
               go to 30
            endif

            wall(1) = 1.d0
            wall(2) = 1.d0
            wall(3) = 1.d0
            if (hR.le.drytol) then
               call riemanntype(hL,hL,uL,-uL,hstar,um,s1m,s2m,&
                                           rare1,rare2,1,drytol,g)
               hstartest=max(hL,hstar)
               ! write(*,*) "slosh from left HIGHER THAN WALL", hstartest+bL.ge.bR
               if (hstartest+bL.lt.bR) then !then !right state should become ghost values that mirror left for wall problem
            !                bR=hstartest+bL
                  wall(2)=0.d0
                  wall(3)=0.d0
                  hR=hL
                  huR=-huL
                  bR=bL
                  phiR=phiL
                  uR=-uL
                  hvR=hvL
                  vR=vL
               elseif (hL+bL.le.bR) then
                  bR=hL+bL
               endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
               call riemanntype(hR,hR,-uR,uR,hstar,um,s1m,s2m,&
                                           rare1,rare2,1,drytol,g)
               hstartest=max(hR,hstar)
               ! write(*,*) "slosh from rght HIGHER THAN WALL", hstartest+bR.ge.bL
               if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
            !               bL=hstartest+bR
                  wall(1)=0.d0
                  wall(2)=0.d0
                  hL=hR
                  huL=-huR
                  bL=bR
                  phiL=phiR
                  uL=-uR
                  hvL=hvR
                  vL=vR
               elseif (hR+bR.le.bL) then
                  bL=hR+bR
               endif
            endif

            !determine wave speeds
            sL=uL-sqrt(g*hL) ! 1 wave speed of left state
            sR=uR+sqrt(g*hR) ! 2 wave speed of right state
            uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
            chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
            sRoe1=uhat-chat ! Roe wave speed 1 wave
            sRoe2=uhat+chat ! Roe wave speed 2 wave

            sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
            sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

            !--------------------end initializing...finally----------
            !solve Riemann problem.

            !maxiter = 1
            fwr = 0.d0
            swr = 0.d0
            pL = 0.d0
            pR = 0.d0

           call riemann_aug_JCP(1,3,3,hL,hR,huL,&
                huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,&
                                            drytol,g,swr,fwr)
      !          call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
      !      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
      !      &     rho,sw,fw)

             !   call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,&
             ! bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,1.d0,swr,&
             !    fwr)

            !eliminate ghost fluxes for wall
            if (swr(2).lt.0.d0) then
              wall(2) = wall(1)
            else if (swr(2).gt.0.d0) then
              wall(2) = wall(3)
            else
              wall(2) = 0.d0
            end if
            ! print *, "WALL:", wall
            do mw=1,3
               s(mw,i)=swr(mw)*wall(mw)
                  fwave(1,mw,i)=fwr(1,mw)*wall(mw)
                  fwave(mu,mw,i)=fwr(2,mw)*wall(mw)
                  fwave(nv,mw,i)=fwr(3,mw)*wall(mw)
            end do
            ! print *, "Spped: ", s(:,i)

 30         continue
            do  mw=1,3
               if (s(mw,i) < 0.d0) then
                     amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
               else if (s(mw,i) > 0.d0) then
                  apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
               else
                 amdq(1:3,i) = amdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
                 apdq(1:3,i) = apdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
               endif
            enddo
          enddo
          ! print *,"AMDQ:", amdq
          s_wall(1) = 0.5d0*(s(1,1)+s(1,2))!
          s_wall(2) = 0.5d0*(s(2,1)+s(2,2))!
          s_wall(3) = 0.5d0*(s(3,1)+s(3,2))!

          if (.not. L2R .and. .not. R2L) then
            s_wall(1) =minval(s)
            s_wall(2) =0.d0
            s_wall(3) =maxval(s)
          end if

            wave_wall = 0.d0

            if ( s_wall(3)-s_wall(1)/=0.d0 ) then
              ! wave_wall(1,1) = (s_wall(3)*s_wall(2)*sum(fwave(1,:,:))&
              ! -(s_wall(2)+s_wall(3))*sum(fwave(mu,:,:)) + sum(fwave(nv,:,:))) /&
              ! ((s_wall(3)-s_wall(1))*(s_wall(2)-s_wall(1)))
              ! wave_wall(1,3) = (-(s_wall(2)+s_wall(1))*sum(fwave(mu,:,:) +&
              ! s_wall(2)*s_wall(1)*sum(fwave(1,:,:)+sum(fwave(nv,:,:))))) /&
              ! ((s_wall(3)-s_wall(1)*(s_wall(3)-s_wall(2))))
              ! wave_wall(1,2) = (s_wall(1)*s_wall(3)*sum(fwave(1,:,:) -&
              ! (s_wall(1)+s_wall(3))*sum(fwave(mu,:,:))+sum(fwave(nv,:,:)))) /&
              ! ((s_wall(2)-s_wall(3)*(s_wall(2)-s_wall(1))))
              ! wave_wall(mu,1) = wave_wall(1,1)*s_wall(1)
              ! wave_wall(nv,1) = wave_wall(2,1)*s_wall(1)
              ! wave_wall(mu,3) = wave_wall(1,3)*s_wall(3)
              ! wave_wall(nv,3) = wave_wall(2,3)*s_wall(3)
              ! wave_wall(mu,2) = wave_wall(1,2)*s_wall(2)
              ! wave_wall(nv,2) = wave_wall(2,2)*s_wall(2)


      wave_wall(1,1)=(s_wall(3)*sum(fwave(1,:,:))-sum(fwave(mu,:,:)))&
        / (s_wall(3)-s_wall(1))
      wave_wall(1,3)=(sum(fwave(mu,:,:))-s_wall(1)*sum(fwave(1,:,:)))&
        / (s_wall(3)-s_wall(1))
      wave_wall(nv,2) = -v_av*sum(fwave(1,:,:)) + sum(fwave(nv,:,:))
            wave_wall(mu,1) = wave_wall(1,1)*s_wall(1)
            wave_wall(nv,1) = wave_wall(1,1)*v_av
            wave_wall(mu,3) = wave_wall(1,3)*s_wall(3)
            wave_wall(nv,3) = wave_wall(1,3)*v_av
            end if

            amdq_wall = 0.d0
            apdq_wall = 0.d0

            do mw=1,size(s_wall)
              if ( s_wall(mw).lt.0.d0 ) then
                amdq_wall(1:3) = amdq_wall(1:3) + wave_wall(:,mw)
              else if (s_wall(mw).gt.0.d0) then
                apdq_wall(1:3) = apdq_wall(1:3) + wave_wall(:,mw)
              else
                amdq_wall(1:3) = amdq_wall(1:3) + 0.5d0*wave_wall(:,mw)
                apdq_wall(1:3) = apdq_wall(1:3) + 0.5d0*wave_wall(:,mw)
              end if
            end do
            !
            ! print *, "AMDQ: ", amdq_wall
            ! print *, "APDQ: ", apdq_wall
            ! print *, "wavw:" ,fwave

          end subroutine redistribute_fwave
    function f(Q,ixy) result(F_vec)
      implicit none
      integer :: ixy,mu,nv
      real(kind=8) :: Q(3),F_vec(3),grav
      common /cparam/  grav
      g = grav
      drytol = dry_tolerance
      if ( Q(1).lt.drytol ) then
        F_vec(1:3) = 0.d0
        return
      else
        if ( ixy==1 ) then
          mu = 2
          nv = 3
        else
          mu = 3
          nv = 2
        end if
        F_vec(1) = Q(mu)
        F_vec(mu) = (Q(mu)**2)/(Q(1)) + 0.5d0*g*Q(1)**2
        F_vec(nv) = Q(nv)*Q(mu)/Q(1)
        return
      end if
    end function f
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function f_rot(q_rot,alpha,beta)
      implicit none
      real(8) :: alpha, beta
      real(8) :: q_rot(3), f_rot(3)
      f_rot = alpha*f(q_rot,1) + beta*f(q_rot,2)
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine riemann_solver(ql,qr,bL,bR,s,fwave,amdq,apdq,ixy)
      implicit none
      integer :: ixy
      real(8) :: ql(3),qr(3),bL,bR,s(3),fwave(3,3),amdq(3),apdq(3)
      ! local
      real(8) :: swl(3),fwl(3,3)
      real(8) :: hL,hR,uL,um,uR,huL,huR,vL,vR,hvL,hvR,wall(3),phiL,phiR
      integer :: mw,mu,nv
      real(8) :: hstar,s1m,s2m,hstartest,sL,sR,uhat,chat,sRoe1,sRoe2
      real(8) :: sE1,sE2,grav,dry_tolerance,g,drytol
      logical :: rare1,rare2
      common /cparam/ grav,dry_tolerance

      g= grav
      drytol = dry_tolerance
      ! print *, g,drytol
      do mw=1,3
         s(mw)=0.d0
        fwave(1,mw)=0.d0
        fwave(2,mw)=0.d0
        fwave(3,mw)=0.d0
      enddo

      ! set normal direction
      if (ixy.eq.1) then
         mu=2
         nv=3
      else
         mu=3
         nv=2
      endif

      ! zero (small) negative values if they exist
      if (qr(1).lt.0.d0) then
            qr(1)=0.d0
            qr(2)=0.d0
            qr(3)=0.d0
      endif

      if (ql(1).lt.0.d0) then
            ql(1)=0.d0
            ql(2)=0.d0
            ql(3)=0.d0
      endif

      !skip problem if in a completely dry area
      if (qr(1) <= drytol .and. ql(1) <= drytol) then
         go to 30
      endif

      !Riemann problem variables
      hL = ql(1)
      hR = qr(1)
      huL = ql(mu)
      huR = qr(mu)
      hvL=ql(nv)
      hvR=qr(nv)

      !check for wet/dry boundary
      if (hR.gt.drytol) then
         uR=huR/hR
         vR=hvR/hR
         phiR = 0.5d0*g*hR**2 + huR**2/hR
      else
         hR = 0.d0
         huR = 0.d0
         hvR = 0.d0
         uR = 0.d0
         vR = 0.d0
         phiR = 0.d0
      endif

      if (hL.gt.drytol) then
         uL=huL/hL
         vL=hvL/hL
         phiL = 0.5d0*g*hL**2 + huL**2/hL
      else
         hL=0.d0
         huL=0.d0
         hvL=0.d0
         uL=0.d0
         vL=0.d0
         phiL = 0.d0
      endif

      wall(1) = 1.d0
      wall(2) = 1.d0
      wall(3) = 1.d0
      if (hR.le.drytol) then
         call riemanntype(hL,hL,uL,-uL,hstar,um,s1m,s2m,rare1,rare2,1,drytol,g)
         hstartest=max(hL,hstar)
         if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
    !                bR=hstartest+bL
            wall(2)=0.d0
            wall(3)=0.d0
            hR=hL
            huR=-huL
            bR=bL
            phiR=phiL
            uR=-uL
            vR=vL
         elseif (hL+bL.lt.bR) then
            bR=hL+bL
         endif
      elseif (hL.le.drytol) then ! right surface is lower than left topo
         call riemanntype(hR,hR,-uR,uR,hstar,um,s1m,s2m,rare1,rare2,1,drytol,g)
         hstartest=max(hR,hstar)
         if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
    !               bL=hstartest+bR
            wall(1)=0.d0
            wall(2)=0.d0
            hL=hR
            huL=-huR
            bL=bR
            phiL=phiR
            uL=-uR
            vL=vR
         elseif (hR+bR.lt.bL) then
            bL=hR+bR
         endif
      endif

      !determine wave speeds
      sL=uL-sqrt(g*hL) ! 1 wave speed of left state
      sR=uR+sqrt(g*hR) ! 2 wave speed of right state

      uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
      chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
      sRoe1=uhat-chat ! Roe wave speed 1 wave
      sRoe2=uhat+chat ! Roe wave speed 2 wave

      sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
      sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

      !--------------------end initializing...finally----------
      !solve Riemann problem.

      call riemann_aug_JCP(1,3,3,hL,hR,huL, &
           huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2, &
                                       drytol,g,swl,fwl)
    !        !eliminate ghost fluxes for wall
      do mw=1,3
         swl(mw)=swl(mw)*wall(mw)
            fwl(1,mw)=fwl(1,mw)*wall(mw)
            fwl(2,mw)=fwl(2,mw)*wall(mw)
            fwl(3,mw)=fwl(3,mw)*wall(mw)
      enddo

      do mw=1,3
         s(mw)=swl(mw)
         fwave(1,mw)=fwl(1,mw)
         fwave(mu,mw)=fwl(2,mw)
         fwave(nv,mw)=fwl(3,mw)
      enddo
    30 continue
      amdq(1:3) = 0.d0
      apdq(1:3) = 0.d0
      do  mw=1,3
         if (s(mw) < 0.d0) then
               amdq(1:3) = amdq(1:3) + fwave(1:3,mw)
         else if (s(mw) > 0.d0) then
            apdq(1:3)  = apdq(1:3) + fwave(1:3,mw)
         else
           amdq(1:3) = amdq(1:3) + 0.5d0 * fwave(1:3,mw)
           apdq(1:3) = apdq(1:3) + 0.5d0 * fwave(1:3,mw)
         endif
      enddo

    end subroutine



    subroutine riemanntype(hL,hR,uL,uR,hm,um,s1m,s2m,rare1,rare2,&
                maxiter,drytol,g)

    !determine the Riemann structure (wave-type in each family)


    implicit none

    !input
    real(8):: hL,hR,uL,uR,drytol,g
    integer maxiter

    !output
    real(8):: s1m,s2m
    logical rare1,rare2

    !local
    real(8):: hm,u1m,u2m,um,delu
    real(8):: h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
    integer iter



    h_min=min(hR,hL)
    h_max=max(hR,hL)
    delu=uR-uL

    if (h_min.le.drytol) then
       hm=0.d0
       um=0.d0
       s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
       s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
       if (hL.le.0.d0) then
          rare2=.true.
          rare1=.false.
       else
          rare1=.true.
          rare2=.false.
       endif

    else
       F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
       F_max= delu + &
            (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

       if (F_min.gt.0.d0) then !2-rarefactions

          hm=(1.d0/(16.d0*g))* &
                  max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
          um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

          s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
          s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

          rare1=.true.
          rare2=.true.

       elseif (F_max.le.0.d0) then !2 shocks

          h0=h_max
          do iter=1,maxiter
             gL=sqrt(.5d0*g*(1/h0 + 1/hL))
             gR=sqrt(.5d0*g*(1/h0 + 1/hR))
             F0=delu+(h0-hL)*gL + (h0-hR)*gR
             dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+ &
                    gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
             slope=2.d0*sqrt(h0)*dfdh
             h0=(sqrt(h0)-F0/slope)**2
          enddo
             hm=h0
             u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
             u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
             um=.5d0*(u1m+u2m)

             s1m=u1m-sqrt(g*hm)
             s2m=u2m+sqrt(g*hm)
             rare1=.false.
             rare2=.false.

       else !one shock one rarefaction
          h0=h_min

          do iter=1,maxiter
             F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
                    + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
             slope=(F_max-F0)/(h_max-h_min)
             h0=h0-F0/slope
          enddo

          hm=h0
          if (hL.gt.hR) then
             um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
             s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
             s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
             rare1=.true.
             rare2=.false.
          else
             s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
             s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
             um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
             rare2=.true.
             rare1=.false.
          endif
       endif
    endif
    return

  end subroutine !riemanntype----------------------------------------------------------------
  subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,&
   hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,&
   sw,fw)

  ! solve shallow water equations given single left and right states
  ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
  ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

  ! To use the original solver call with maxiter=1.

  ! This solver allows iteration when maxiter > 1. The iteration seems to help with
  ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
  ! due to loss of hyperbolicity.

  implicit none

  !input
  integer meqn,mwaves,maxiter
  real(8):: sw(3),fw(3,3)
  real(8):: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
  real(8):: hvL,hvR,vL,vR,pL,pR,um
  real(8):: drytol,g,rho


  !local
  integer m,mw,k,iter
  real(8):: A(3,3)
  real(8):: r(3,3)
  real(8):: lambda(3)
  real(8):: del(3)
  real(8):: beta(3)

  real(8):: delh,delhu,delphi,delb,delnorm
  real(8):: rare1st,rare2st,sdelta,raremin,raremax
  real(8):: criticaltol,convergencetol,raretol
  real(8):: criticaltol_2, hustar_interface
  real(8):: s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
  real(8):: huRstar,huLstar,uRstar,uLstar,hstarHLL
  real(8):: deldelh,deldelphi,delP
  real(8):: s1m,s2m,hm
  real(8):: det1,det2,det3,determinant

  logical rare1,rare2,rarecorrector,rarecorrectortest,sonic
  !determine del vectors
  delh = hR-hL
  delhu = huR-huL
  delphi = phiR-phiL
  delb = bR-bL
  delnorm = delh**2 + delphi**2

  call riemanntype(hL,hR,uL,uR,hm,um,s1m,s2m,rare1,rare2,&
                                          5,drytol,g)
  ! print *, "I make it here tho",sE1,s2m,sE2,s1m

  lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
  lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
  sE1=lambda(1)
  sE2=lambda(3)
  lambda(2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

  ! print *, "lambda:", lambda

  hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

!     !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.false.
      rarecorrector=.false.
      if (rarecorrectortest) then
         sdelta=lambda(3)-lambda(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
            rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and. &
              max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  lambda(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  lambda(2)=s2m
               else
                  lambda(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

      do mw=1,mwaves
         r(1,mw)=1.d0
         r(2,mw)=lambda(mw)
         r(3,mw)=(lambda(mw))**2
      enddo
      if (.not.rarecorrector) then
         lambda(2) = 0.5d0*(lambda(1)+lambda(3))
         r(1,2)=0.d0
         r(2,2)=0.d0
         r(3,2)=1.d0
      endif


      !criticaltol = 1.d-6
      ! MODIFIED:
      criticaltol = max(drytol*g, 1d-6)
      criticaltol_2 = sqrt(criticaltol)
      deldelh = -delb
      deldelphi = -0.5d0 * (hR + hL) * (g * delb ) !+ delp / rho)

      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,maxiter
         !determine steady state wave (this will be subtracted from the delta vectors)
         if (min(hLstar,hRstar).lt.drytol.and.rarecorrector) then
            rarecorrector=.false.
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            lambda(2) = 0.5d0*(lambda(1)+lambda(3))
            r(1,2)=0.d0
            r(2,2)=0.d0
            r(3,2)=1.d0
         endif

         hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
         s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
         s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

         ! MODIFIED from 5.3.1 version
         sonic = .false.
         if (abs(s1s2bar) <= criticaltol) then
            sonic = .true.
         else if (s1s2bar*s1s2tilde <= criticaltol**2) then
            sonic = .true.
         else if (s1s2bar*sE1*sE2 <= criticaltol**2) then
            sonic = .true.
         else if (min(abs(sE1),abs(sE2)) < criticaltol_2) then
            sonic = .true.
         else if (sE1 <  criticaltol_2 .and. s1m > -criticaltol_2) then
            sonic = .true.
         else if (sE2 > -criticaltol_2 .and. s2m <  criticaltol_2) then
            sonic = .true.
         else if ((uL+dsqrt(g*hL))*(uR+dsqrt(g*hR)) < 0.d0) then
            sonic = .true.
         else if ((uL- dsqrt(g*hL))*(uR- dsqrt(g*hR)) < 0.d0) then
            sonic = .true.
         end if

         if (sonic) then
            deldelh =  -delb
         else
            deldelh = delb*g*hbar/s1s2bar
         endif
         if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
         elseif (sE1.ge.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
            deldelh = max(deldelh,-hL)
         elseif (sE2.le.-criticaltol) then
            deldelh = min(deldelh,hR)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
         endif

         deldelh = deldelh !- delP/(rho*g)

         if (sonic) then
            deldelphi = -g*hbar*delb
         else
            deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
         endif
         deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
         deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))
         ! deldelphi = deldelphi - hbar * delp / rho

         del(1)=delh-deldelh
         del(2)=delhu
         del(3)=delphi-deldelphi

         det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
         det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
         det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
         determinant=det1-det2+det3

         do k=1,3
            do mw=1,3
                  A(1,mw)=r(1,mw)
                  A(2,mw)=r(2,mw)
                  A(3,mw)=r(3,mw)
            enddo
            A(1,k)=del(1)
            A(2,k)=del(2)
            A(3,k)=del(3)
            det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            beta(k)=(det1-det2+det3)/determinant
         enddo

         !exit if things aren't changing
         if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
         delnorm = del(1)**2+del(3)**2
         !find new states qLstar and qRstar on either side of interface
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar
         do mw=1,mwaves
            if (lambda(mw).lt.0.d0) then
               hLstar= hLstar + beta(mw)*r(1,mw)
               huLstar= huLstar + beta(mw)*r(2,mw)
            endif
         enddo
         do mw=mwaves,1,-1
            if (lambda(mw).gt.0.d0) then
               hRstar= hRstar - beta(mw)*r(1,mw)
               huRstar= huRstar - beta(mw)*r(2,mw)
            endif
         enddo

         if (hLstar.gt.drytol) then
            uLstar=huLstar/hLstar
         else
            hLstar=max(hLstar,0.d0)
            uLstar=0.d0
         endif
         if (hRstar.gt.drytol) then
            uRstar=huRstar/hRstar
         else
            hRstar=max(hRstar,0.d0)
            uRstar=0.d0
         endif

      enddo ! end iteration on Riemann problem
      do mw=1,mwaves
         sw(mw)=lambda(mw)
         fw(1,mw)=beta(mw)*r(2,mw)
         fw(2,mw)=beta(mw)*r(3,mw)
         fw(3,mw)=beta(mw)*r(2,mw)
      enddo

      !find transverse components (ie huv jumps).
      ! MODIFIED from 5.3.1 version
      fw(3,1)=fw(3,1)*vL
      fw(3,3)=fw(3,3)*vR
      fw(3,2)=hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3)

      ! print *, "Fw:",fw

      return

    end subroutine !riemann_aug_JCP-------------------------------------------------
      subroutine riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,&
                 bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,s1,s2,drytol,g,rho,&
                 sw,fw)

      ! solve shallow water equations given single left and right states
      ! solution has two waves.
      ! flux - source is decomposed.

      implicit none

      !input
      integer meqn,mwaves

      real(8):: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,s1,s2
      real(8):: hvL,hvR,vL,vR,pL,pR
      real(8):: drytol,g,rho

      real(8):: sw(mwaves)
      real(8):: fw(meqn,mwaves)

      !local
      real(8):: delh,delhu,delphi,delb,delhdecomp,delphidecomp
      real(8):: deldelh,deldelphi,delP
      real(8):: beta1,beta2


      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL
      ! delP = pR - pL

      deldelphi = -0.5d0 * (hR + hL) * (g * delb) ! + delP / rho)
      delphidecomp = delphi - deldelphi

      !flux decomposition
      beta1 = (s2*delhu - delphidecomp)/(s2-s1)
      beta2 = (delphidecomp - s1*delhu)/(s2-s1)

      sw(1)=s1
      sw(2)=0.5d0*(s1+s2)
      sw(3)=s2
      ! 1st nonlinear wave
      fw(1,1) = beta1
      fw(2,1) = beta1*s1
      fw(3,1) = beta1*vL
      ! 2nd nonlinear wave
      fw(1,3) = beta2
      fw(2,3) = beta2*s2
      fw(3,3) = beta2*vR
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL -fw(3,1)-fw(3,3)
      return

      end subroutine ! -------------------------------------------------
      ! =====================================================
      subroutine rpt22(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      ! =====================================================

            implicit none
      !
      !     # Riemann solver in the transverse direction using an einfeldt
      !     Jacobian.

            double precision :: grav, g
            double precision, parameter :: tol = 1.e-14
            common /cparam/ grav

            integer ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

            double precision  ql(meqn,1-mbc:maxm+mbc)
            double precision  qr(meqn,1-mbc:maxm+mbc)
            double precision  asdq(meqn,1-mbc:maxm+mbc)
            double precision  bmasdq(meqn,1-mbc:maxm+mbc)
            double precision  bpasdq(meqn,1-mbc:maxm+mbc)
            double precision  aux1(maux,1-mbc:maxm+mbc)
            double precision  aux2(maux,1-mbc:maxm+mbc)
            double precision  aux3(maux,1-mbc:maxm+mbc)

            double precision  s(3)
            double precision  r(3,3)
            double precision  beta(3)
            double precision  abs_tol
            double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
            double precision  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
            double precision  delf1,delf2,delf3,dxdcd,dxdcu
            double precision  dxdcm,dxdcp,topo1,topo3,eta

            integer i,m,mw,mu,mv

            abs_tol=tol
            g = grav

            if (ixy.eq.1) then
              mu = 2
              mv = 3
            else
              mu = 3
              mv = 2
            endif

            do i=2-mbc,mx+mbc

               hl=qr(1,i-1)
               hr=ql(1,i)
               hul=qr(mu,i-1)
               hur=ql(mu,i)
               hvl=qr(mv,i-1)
               hvr=ql(mv,i)

      !===========determine velocity from momentum===========================
             if (hl.lt.abs_tol) then
                hl=0.d0
                ul=0.d0
                vl=0.d0
             else
                ul=hul/hl
                vl=hvl/hl
             endif

             if (hr.lt.abs_tol) then
                hr=0.d0
                ur=0.d0
                vr=0.d0
             else
                ur=hur/hr
                vr=hvr/hr
             endif

             do mw=1,mwaves
                s(mw)=0.d0
                beta(mw)=0.d0
                do m=1,meqn
                   r(m,mw)=0.d0
                enddo
             enddo
            dxdcp = 1.d0
            dxdcm = 1.d0

            if (hl <= tol .and. hr <= tol) go to 90

             !check and see if cell that transverse waves are going in is high and dry
             if (imp.eq.1) then
                  eta = qr(1,i-1)  + aux2(1,i-1)
                  topo1 = aux1(1,i-1)
                  topo3 = aux3(1,i-1)
      !            s1 = vl-sqrt(g*hl)
      !            s3 = vl+sqrt(g*hl)
      !            s2 = 0.5d0*(s1+s3)
             else
                  eta = ql(1,i) + aux2(1,i)
                  topo1 = aux1(1,i)
                  topo3 = aux3(1,i)
      !            s1 = vr-sqrt(g*hr)
      !            s3 = vr+sqrt(g*hr)
      !            s2 = 0.5d0*(s1+s3)
             endif
             if (eta.lt.max(topo1,topo3)) go to 90


      !=====Determine some speeds necessary for the Jacobian=================
                  vhat=(vr*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + &
                    (vl*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))

                  uhat=(ur*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + &
                    (ul*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
                  hhat=(hr+hl)/2.d0

                  roe1=vhat-dsqrt(g*hhat)
                  roe3=vhat+dsqrt(g*hhat)

                  s1l=vl-dsqrt(g*hl)
                  s3r=vr+dsqrt(g*hr)

                  s1=dmin1(roe1,s1l)
                  s3=dmax1(roe3,s3r)

                  s2=0.5d0*(s1+s3)

                 s(1)=s1
                 s(2)=s2
                 s(3)=s3
      !=======================Determine asdq decomposition (beta)============
               delf1=asdq(1,i)
               delf2=asdq(mu,i)
               delf3=asdq(mv, i)

               beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
               beta(2) = -s2*delf1 + delf2
               beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
      !======================End =================================================

      !=====================Set-up eigenvectors===================================
               r(1,1) = 1.d0
               r(2,1) = s2
               r(3,1) = s1

               r(1,2) = 0.d0
               r(2,2) = 1.d0
               r(3,2) = 0.d0

               r(1,3) = 1.d0
               r(2,3) = s2
               r(3,3) = s3
      !============================================================================
      90      continue
      !============= compute fluctuations==========================================

                     bmasdq(1,i)=0.0d0
                     bpasdq(1,i)=0.0d0
                     bmasdq(2,i)=0.0d0
                     bpasdq(2,i)=0.0d0
                     bmasdq(3,i)=0.0d0
                     bpasdq(3,i)=0.0d0
                  do  mw=1,3
                     if (s(mw).lt.0.d0) then
                       bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                       bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                       bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
                     elseif (s(mw).gt.0.d0) then
                       bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                       bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                       bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
                     endif
                  enddo
      !========================================================================
               enddo
      !

      !

            return
            end


end module redistribute2D

!     ==========================================================
    subroutine step2ds(maxm,num_eqn,num_waves,num_aux,num_ghost,mx,my, &
                        qold,qnew,qold2,qnew2,aux,auxu,dx,dy,dt,method,&
                        mthlim,cfl,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d, &
                        aux1,aux2l,aux3,work,mwork,ids,use_fwave,rpn2,rpt2,bar_ht)
!     ==========================================================

!     # Take one time step, updating q.
!     # On entry, qold and qnew should be identical and give the
!     #    initial data for this step
!     # On exit, qnew returns values at the end of the time step.
!     #    qold is unchanged.

!     # qadd is used to return increments to q from flux2
!     # fadd and gadd are used to return flux increments from flux2.
!     # See the flux2 documentation for more information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use redistribute2D
    use aux_module_SRD
    implicit double precision (a-h,o-z)
    double precision :: qold(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qold2(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qnew(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qnew2(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision ::  q1d(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision ::  q1d2(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: qadd(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: qadd2(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: fadd(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: fadd2(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: gadd(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
    double precision :: gadd2(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
    dimension bmasdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension bpasdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension bmasdq2(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension bpasdq2(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: auxu(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: aux1(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux2l(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux2u(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux3(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: dtdx1d(1-num_ghost:maxm+num_ghost)
    double precision :: dtdy1d(1-num_ghost:maxm+num_ghost)
    integer :: method(7),mthlim(num_waves)
    logical ::          use_fwave
    double precision :: work(mwork)
    double precision :: q_final(3,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension fp1(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension fm1(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension gp1(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension gm1(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension fp2(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension fm2(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension gp2(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension gm2(num_eqn,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    dimension apdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension amdq(num_eqn, 1-num_ghost:maxm+num_ghost)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! barrier parameters
    integer :: N, ii(mx*4),jj(mx*4),type_supper(mx*4),type_sunder(mx*4)
    integer :: N2, ii2(mx*4),jj2(mx*4),type_supper2(mx*4),type_sunder2(mx*4)
    integer :: k_count_up,all_undercells_i(mx*4),all_undercells_j(mx*4)
    integer :: k_count_un,N_ij_un(1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    integer :: N_ij_up(1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    integer :: all_uppercells_i(mx*4)
    integer :: all_uppercells_j(mx*4),unS_cells_i(mx*4),unS_cells_j(mx*4)
    integer :: upS_cells_i(mx*4),upS_cells_j(mx*4),unS_cells_i2(mx*4),unS_cells_j2(mx*4)
    integer :: upS_cells_i2(mx*4),upS_cells_j2(mx*4)
    real(8) :: up_area_ij(1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    real(8) :: un_area_ij(1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost)
    real(8) :: intersections(2,4*mx),lengths_sunder(5,mx*4),lengths_supper(5,mx*4)
    real(8) :: intersections2(2,4*mx),lengths_sunder2(5,mx*4),lengths_supper2(5,mx*4)
    real(8) :: area_supper(mx*4), area_sunder(mx*4),xlower,ylower,xupper,yupper
    real(kind=8)::x_0,x_e,y_0,y_e,alpha,bar_height,bar_ht,area_supper2(mx*4), area_sunder2(mx*4)
    integer :: i_0,i_e,j_0,j_e,bar_index,bar_index_i,k
    integer :: ind(mx+2),ind2(mx+2), iostatus,kk
    logical :: L2R,R2L,xor_lr,check_on,lexist,ot
    real(kind=8) :: xe(1-num_ghost:mx+num_ghost+1), ye(1-num_ghost:my+num_ghost+1)
    integer :: i1,j1,i2,j2
    real :: a1,a2,cfl2,qval(3),qval2(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external :: rpn2,rpt2
    ! barrier param set on setup python file
    ! common /cparam/ bar_index_i,bar_ht,bar_loc
!f2py intent(out) cfl
!f2py intent(in,out) qnew, qnew2
!f2py optional q1d, qadd, fadd, gadd, dtdx1d, dtdy1d

! Dummy interfaces just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn2(x)
!f2py x=rpt2(x)

!     # partition work array into pieces needed for local storage in
!     # flux2 routine.  Find starting index of each piece:

    i0wave = 1
    i0s = i0wave + (maxm+2*num_ghost)*num_eqn*num_waves
    i0amdq = i0s + (maxm+2*num_ghost)*num_waves
    i0apdq = i0amdq + (maxm+2*num_ghost)*num_eqn
    i0cqxx = i0apdq + (maxm+2*num_ghost)*num_eqn
    i0bmadq = i0cqxx + (maxm+2*num_ghost)*num_eqn
    i0bpadq = i0bmadq + (maxm+2*num_ghost)*num_eqn
    iused = i0bpadq + (maxm+2*num_ghost)*num_eqn - 1

    if (iused > mwork) then
    !        # This shouldn't happen due to checks in claw2
        write(6,*) '*** not enough work space in step2'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif
    ! get the preprocessed data of the geometry of small cells:
    open (unit=1,file="./small_cells_data.txt",FORM="FORMATTED",STATUS="OLD",&
    ACTION="READ",access='sequential')
    rewind 1
    read (1,*,end=101) N
    read (1,*) N2
    read (1,*) ii(1:N)
    read (1,*) jj(1:N)
    read (1,*) ii2(1:N2)
    read (1,*) jj2(1:N2)
    read (1,*) intersections(:,1:N+1)
    read (1,*) intersections2(:,1:N2+1)
    read (1,*) type_supper(1:N)
    read (1,*) type_sunder(1:N)
    read (1,*) type_supper2(1:N2)
    read (1,*) type_sunder2(1:N2)
    read (1,*) lengths_supper(:,1:N)
    read (1,*) lengths_sunder(:,1:N)
    read (1,*) lengths_supper2(:,1:N2)
    read (1,*) lengths_sunder2(:,1:N2)
    read (1,*) unS_cells_i(1:N)
    read (1,*) unS_cells_j(1:N)
    read (1,*) unS_cells_i2(1:N2)
    read (1,*) unS_cells_j2(1:N2)
    read (1,*) upS_cells_i(1:N)
    read (1,*) upS_cells_j(1:N)
    read (1,*) upS_cells_i2(1:N2)
    read (1,*) upS_cells_j2(1:N2)
    read (1,*) area_supper(1:N)
    read (1,*) area_sunder(1:N)
    read (1,*) area_supper2(1:N2)
    read (1,*) area_sunder2(1:N2)
    read (1,*) xlower
    read (1,*) ylower
    read (1,*) xupper
    read (1,*) yupper
    read (1,*) x_0
    read (1,*) y_0
    read (1,*) x_1
    read (1,*) y_1
    read (1,*) x_e
    read (1,*) y_e
    read (1,*) x_2
    read (1,*) y_2
    read (1,*) up_area_ij
    read (1,*) un_area_ij
    ! ENDFILE 1
    ! REWIND 1
    close (1,status="keep")

  101 continue
    index_capa = method(6)
    num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! grid and barrier info
    bar_height = bar_ht
    ! print*,"bar ht:",bar_ht,bar_height
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (index_capa == 0) then
    !        # no capa array:
        do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
        5 END DO
    endif
    fm1=0.d0
    fp1=0.d0
    fm2=0.d0
    fp2=0.d0
    gm1 =0.d0
    gp1=0.d0
    gm2=0.d0
    gp2=0.d0
    forall (m=1:num_eqn, i = 1-num_ghost: mx+num_ghost)
    gadd(m,1,i) = 0.d0
    gadd(m,2,i) = 0.d0
    gadd2(m,1,i) = 0.d0
    gadd2(m,2,i) = 0.d0
    end forall
    ! if( ids == 1 )then

    !     # perform x-sweeps
    !     ==================

    !     # note that for dimensional splitting we sweep over the rows of
    !     # ghosts cells as well as the interior.  This updates the ghost
    !     # cell values to the intermediate state as needed in the following
    !     # sweep in the y-direction.

        do 50 j = 0,my+1

        !        # copy data along a slice into 1d arrays:
            forall (m=1:num_eqn, i = 1-num_ghost: mx+num_ghost)
            q1d(m,i) = qold(m,i,j)
            q1d2(m,i) = qold2(m,i,j)
            end forall
            ! print *, "Q1D XSLICE:", q1d(1,:)
            ! print*, "Q1D2 XSLICE:", q1d2(1,:)
            if (index_capa > 0)  then
                do 22 i = 1-num_ghost, mx+num_ghost
                    dtdx1d(i) = dtdx / aux(index_capa,i,j)
                22 END DO
            endif

            if (num_aux > 0)  then
                do 23 ma=1,num_aux
                    do 23 i = 1-num_ghost, mx+num_ghost
                        aux2l(ma,i) = aux(ma,i,j  )
                        aux2u(ma,i) = auxu(ma,i,j  )
                23 END DO

                if(j /= 1-num_ghost)then
                    do 24 ma=1,num_aux
                        do 24 i = 1-num_ghost, mx+num_ghost
                            aux1(ma,i) = aux(ma,i,j-1)
                    24 END DO
                endif

                if(j /= my+num_ghost)then
                    do 25 ma=1,num_aux
                        do 25 i = 1-num_ghost, mx+num_ghost
                            aux3(ma,i) = aux(ma,i,j+1)
                    25 END DO
                endif

            endif

        !        # Store the value of j along this slice in the common block
        !        # comxyt in case it is needed in the Riemann solver (for
        !        # variable coefficient problems)
            jcom = j

        !        # compute modifications fadd and gadd to fluxes along this slice:
            call flux2(1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx, &
            q1d,dtdx1d,aux1,aux2l,aux3,method,mthlim, &
            qadd,fadd,gadd,cfl1d, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave)
            fp1(:,:,j) = apdq
            fm1(:,:,j) = amdq
            ! print*, "q1d2: ",q1d2
            call rpt22(1,1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d,q1d, &
            aux1,aux2l,aux3,amdq,bmasdq,bpasdq)

            !     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
            forall ( m=1:num_eqn , i = 1:mx+1)
            gadd(m,1,i-1) = gadd(m,1,i-1) - &
            0.5d0*dtdx1d(i-1) * bmasdq(m,i)
            gadd(m,2,i-1) = gadd(m,2,i-1) - &
            0.5d0*dtdx1d(i-1) * bpasdq(m,i)
            end forall
            call rpt22(1,2,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d,q1d, &
            aux1,aux2l,aux3,apdq,bmasdq,bpasdq)

            !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
            forall ( m=1:num_eqn , i = 1: mx+1)
            gadd(m,1,i) = gadd(m,1,i) - &
            0.5d0*dtdx1d(i) * bmasdq(m,i)
            gadd(m,2,i) = gadd(m,2,i) - &
            0.5d0*dtdx1d(i) * bpasdq(m,i)
            end forall


            call flux2(1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx, &
            q1d2,dtdx1d,aux1,aux2u,aux3,method,mthlim, &
            qadd2,fadd2,gadd2,cfl1d2, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave)
            fp2(:,:,j) = apdq
            fm2(:,:,j) = amdq
            ! print*, "apdq",apdq,"amdq",amdq
            cfl = dmax1(cfl,min(cfl1d,cfl1d2))
            ! cfl = dmax1(cfl,cfl1d2)
        !     # split the right-going flux difference into down-going and up-going:
            call rpt22(1,1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d2,q1d2, &
            aux1,aux2u,aux3,amdq,bmasdq2,bpasdq2)

        !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
        forall ( m=1:num_eqn , i = 1:mx+1)
        gadd2(m,1,i-1) = gadd2(m,1,i-1) - &
        0.5d0*dtdx1d(i-1) * bmasdq2(m,i)
        gadd2(m,2,i-1) = gadd2(m,2,i-1) - &
        0.5d0*dtdx1d(i-1) * bpasdq2(m,i)
        end forall
            call rpt22(1,2,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d2,q1d2, &
            aux1,aux2u,aux3,apdq,bmasdq2,bpasdq2)

        !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
            forall ( m=1:num_eqn , i = 1: mx+1)
            gadd2(m,1,i) = gadd2(m,1,i) - &
            0.5d0*dtdx1d(i) * bmasdq2(m,i)
            gadd2(m,2,i) = gadd2(m,2,i) - &
            0.5d0*dtdx1d(i) * bpasdq2(m,i)
            end forall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! if( ids == 1 )then
            if (index_capa == 0) then

            !            # no capa array.  Standard flux differencing:
                ! forall (m=1:num_eqn, i=1:mx)
                do m=1,num_eqn
                  do i=1,mx
                    qnew(m,i,j) = qnew(m,i,j) + (qadd(m,i) &
                    - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                    - dtdx * (gadd(m,2,i) - gadd(m,1,i)))
                    do i1=1,N
                      if (ii(i1).eq.i .and. jj(i1).eq.j) then
                        goto 1111
                      end if
                    end do
                    do i1=1,N2
                      if (ii2(i1).eq.i .and. jj2(i1).eq.j) then
                        goto 1111
                      end if
                    end do
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i)

                    1111 continue


                    qnew2(m,i,j) = qnew2(m,i,j) + (qadd2(m,i) &
                    - dtdx * (fadd2(m,i+1) - fadd2(m,i)) &
                    - dtdx * (gadd2(m,2,i) - gadd2(m,1,i)))
                    do i1=1,N
                      if (ii(i1).eq.i .and. jj(i1).eq.j) then
                        goto 111
                      end if
                    end do
                    do i1=1,N2
                      if (ii2(i1).eq.i .and. jj2(i1).eq.j) then
                        goto 111
                      end if
                    end do
                    qnew2(m,i,j+1) = qnew2(m,i,j+1) + dtdy * gadd2(m,2,i)
                    qnew2(m,i,j-1) = qnew2(m,i,j-1) - dtdy * gadd2(m,1,i)

                        111        continue


                      ! qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                      ! - dtdx * (fadd(m,i+1) - fadd(m,i))
                      ! qnew2(m,i,j) = qnew2(m,i,j) + qadd2(m,i) &
                      ! - dtdx * (fadd2(m,i+1) - fadd2(m,i))
                  end do
                end do
                ! end forall
            else
              print *, "I here?"

            !            # with capa array.
                forall (m=1:num_eqn, i=1:mx)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                / auxu(index_capa,i,j)
                end forall
            endif
            ! endif
        50 END DO
    ! endif
    !
    ! if( ids == 2 )then
    forall (m=1:num_eqn, i = 1-num_ghost: mx+num_ghost)
    gadd(m,1,i) = 0.d0
    gadd(m,2,i) = 0.d0
    gadd2(m,1,i) = 0.d0
    gadd2(m,2,i) = 0.d0
    end forall
    !     # perform y sweeps
    !     ==================

        do 100 i = 0, mx+1

        !        # copy data along a slice into 1d arrays:
            forall (m=1:num_eqn, j = 1-num_ghost: my+num_ghost)
            q1d(m,j) = qold(m,i,j)
            q1d2(m,j) = qold2(m,i,j)
            end forall
            ! print *, "Q1D YSLICE:", q1d(1,:)
            ! print*, "Q1D2 YSLICE:", q1d2(1,:)
            if (index_capa > 0)  then
              print *, "I here?"

                do 72 j = 1-num_ghost, my+num_ghost
                    dtdy1d(j) = dtdy / aux(index_capa,i,j)
                72 END DO
            endif

            if (num_aux > 0)  then

                do 73 ma=1,num_aux
                    do 73 j = 1-num_ghost, my+num_ghost
                        aux2l(ma,j) = aux(ma,i,j)
                        aux2u(ma,j) = auxu(ma,i,j)
                73 END DO

                if(i /= 1-num_ghost)then
                    do 74 ma=1,num_aux
                        do 74 j = 1-num_ghost, my+num_ghost
                            aux1(ma,j) = aux(ma,i-1,j)
                    74 END DO
                endif

                if(i /= mx+num_ghost)then
                    do 75 ma=1,num_aux
                        do 75 j = 1-num_ghost, my+num_ghost
                            aux3(ma,j) = aux(ma,i+1,j)
                    75 END DO
                endif

            endif

        !     # Store the value of i along this slice in the common block
        !        # comxyt in case it is needed in the Riemann solver (for
        !        # variable coefficient problems)
            icom = i

        !        # compute modifications fadd and gadd to fluxes along this slice:
            call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my, &
            q1d,dtdy1d,aux1,aux2l,aux3,method,mthlim, &
            qadd,fadd,gadd,cfl1d, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave)
            gm1(:,i,:) = amdq
            gp1(:,i,:) = apdq

            call rpt22(2,1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d,q1d, &
            aux1,aux2l,aux3,amdq,bmasdq,bpasdq)

            !     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
            forall ( m=1:num_eqn , i = 1:mx+1)
            gadd(m,1,i-1) = gadd(m,1,i-1) - &
            0.5d0*dtdx1d(i-1) * bmasdq(m,i)
            gadd(m,2,i-1) = gadd(m,2,i-1) - &
            0.5d0*dtdx1d(i-1) * bpasdq(m,i)
            end forall
            call rpt22(2,2,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d,q1d, &
            aux1,aux2l,aux3,apdq,bmasdq,bpasdq)

            !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
            forall ( m=1:num_eqn , i = 1: mx+1)
            gadd(m,1,i) = gadd(m,1,i) - &
            0.5d0*dtdx1d(i) * bmasdq(m,i)
            gadd(m,2,i) = gadd(m,2,i) - &
            0.5d0*dtdx1d(i) * bpasdq(m,i)
            end forall



            call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my, &
            q1d2,dtdy1d,aux1,aux2u,aux3,method,mthlim, &
            qadd2,fadd2,gadd2,cfl1d2, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave)
            gm2(:,i,:) = amdq
            gp2(:,i,:) = apdq

            cfl = dmax1(cfl,min(cfl1d,cfl1d2))
            ! cfl =dmax1(cfl,cfl1d2)

            call rpt22(2,1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d2,q1d2, &
            aux1,aux2u,aux3,amdq,bmasdq2,bpasdq2)

        !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
            forall ( m=1:num_eqn , i = 1:mx+1)
            gadd2(m,1,i-1) = gadd2(m,1,i-1) - &
            0.5d0*dtdx1d(i-1) * bmasdq2(m,i)
            gadd2(m,2,i-1) = gadd2(m,2,i-1) - &
            0.5d0*dtdx1d(i-1) * bpasdq2(m,i)
            end forall
            call rpt22(2,2,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,q1d2,q1d2, &
            aux1,aux2u,aux3,apdq,bmasdq2,bpasdq2)

            !     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
            forall ( m=1:num_eqn , i = 1: mx+1)
            gadd2(m,1,i) = gadd2(m,1,i) - &
            0.5d0*dtdx1d(i) * bmasdq2(m,i)
            gadd2(m,2,i) = gadd2(m,2,i) - &
            0.5d0*dtdx1d(i) * bpasdq2(m,i)
            end forall
        !        # update qnew by flux differencing.
        !        # Note that the roles of fadd and gadd are reversed for
        !        # the y-sweeps -- fadd is the modification to g-fluxes and
        !        # gadd is the modification to f-fluxes to the left and right.
      ! if( ids == 2 )then
            if (index_capa == 0) then
              ! forall (m=1:num_eqn,j=1:my)
                do m=1,num_eqn
                  do j=1,my

                    qnew(m,i,j) = qnew(m,i,j) + (qadd(m,j) &
                    - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                    - dtdx * (gadd(m,2,j) - gadd(m,1,j)))
                    do i1=1,N
                      if (ii(i1).eq.i .and. jj(i1).eq.j) then
                        goto 1122
                      end if
                    end do
                    do i1=1,N2
                      if (ii2(i1).eq.i .and. jj2(i1).eq.j) then
                        goto 1122
                      end if
                    end do
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j)

                        1122 continue

                    qnew2(m,i,j) = qnew2(m,i,j) + (qadd2(m,j) &
                    - dtdy * (fadd2(m,j+1) - fadd2(m,j)) &
                    - dtdx * (gadd2(m,2,j) - gadd2(m,1,j)))
                    do i1=1,N
                      if (ii(i1).eq.i .and. jj(i1).eq.j) then
                        goto 112
                      end if
                    end do
                    do i1=1,N2
                      if (ii2(i1).eq.i .and. jj2(i1).eq.j) then
                        goto 112
                      end if
                    end do
                    qnew2(m,i+1,j) = qnew2(m,i+1,j) + dtdx * gadd2(m,2,j)
                    qnew2(m,i-1,j) = qnew2(m,i-1,j) - dtdx * gadd2(m,1,j)

                        112          continue



                      ! qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                      ! - dtdx * (fadd(m,i+1) - fadd(m,i))
                      ! qnew2(m,i,j) = qnew2(m,i,j) + qadd2(m,j) &
                      ! - dtdx * (fadd2(m,i+1) - fadd2(m,i))
                  end do
                end do
                ! end forall

            else
              print *, "I here?"

            !            # with capa array.
                forall (m=1:num_eqn, j=1:my)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                / auxu(index_capa,i,j)
                end forall

            endif

        ! endif

        100 END DO


    ! endif

    ! Correct to get the Q^ for both sides of the barrier, for small cells
      call small_cells_update_correct(ii,jj,dtdx,qold,qold2,qnew,qnew2,aux(1,:,:),&
      intersections(:,1:N+1),type_supper(1:N),type_sunder(1:N),lengths_supper(:,1:N),&
      lengths_sunder(:,1:N),area_sunder(1:N),area_supper(1:N),1,N,xlower,ylower,&
      xupper,yupper,mx,my,num_ghost,x_0,y_0,x_e,y_e,bar_height,fm1,fp1,gm1,gp1,fm2,&
      fp2,gm2,gp2,ot)
      call small_cells_update_correct(ii2,jj2,dtdx,qold,qold2,qnew,qnew2,aux(1,:,:),&
      intersections2(:,1:N2+1),type_supper2(1:N2),type_sunder2(1:N2),lengths_supper2(:,1:N2),&
      lengths_sunder2(:,1:N2),area_sunder2(1:N2),area_supper2(1:N2),2,N2,xlower,ylower,&
      xupper,yupper,mx,my,num_ghost,x_1,y_1,x_2,y_2,bar_height,fm1,fp1,gm1,gp1,fm2,&
      fp2,gm2,gp2,ot)
      ! if (ids.eq.2) then

      call SRD_update_correct(qnew,qnew2,all_undercells_i(1:k_count_un),&
                      all_undercells_j(1:k_count_un),all_uppercells_i(1:k_count_up),&
                      all_uppercells_j(1:k_count_up),ii(1:N),jj(1:N),N_ij_up,N_ij_un,&
                      unS_cells_i(1:N),unS_cells_j(1:N),upS_cells_i(1:N),upS_cells_j(1:N),&
                      up_area_ij,un_area_ij,mx,my,num_ghost)
      call SRD_update_correct(qnew,qnew2,all_undercells_i(1:k_count_un),&
                      all_undercells_j(1:k_count_un),all_uppercells_i(1:k_count_up),&
                      all_uppercells_j(1:k_count_up),ii2(1:N2),jj2(1:N2),N_ij_up,N_ij_un,&
                      unS_cells_i2(1:N2),unS_cells_j2(1:N2),upS_cells_i2(1:N2),upS_cells_j2(1:N2),&
                      up_area_ij,un_area_ij,mx,my,num_ghost)
    ! if (ids == 2) then
    !   ! call small_cells_update_correct(ii,jj,dtdx,qold,qold2,qnew,qnew2,aux(1,:,:),&
    !   !   intersections(:,1:N+1),type_supper(1:N),type_sunder(1:N),lengths_supper(:,1:N),&
    !   !   lengths_sunder(:,1:N),area_sunder(1:N),area_supper(1:N),ids,N,xlower,ylower,&
    !   !   xupper,yupper,mx,my,num_ghost,x_0,y_0,x_e,y_e,bar_height,fm1,fp1,gm1,gp1,fm2,&
    !   !   fp2,gm2,gp2)
    ! call SRD_update_correct(qnew,qnew2,all_undercells_i(1:k_count_un),&
    !                 all_undercells_j(1:k_count_un),all_uppercells_i(1:k_count_up),&
    !                 all_uppercells_j(1:k_count_up),ii(1:N),jj(1:N),N_ij_up,N_ij_un,&
    !                 unS_cells_i(1:N),unS_cells_j(1:N),upS_cells_i(1:N),upS_cells_j(1:N),&
    !                 up_area_ij,un_area_ij,mx,my,num_ghost)
    ! ! Correct to get the SRD update using neighborhood values:

    !   ! do i=1,mx
    !   !   do j=1,my
    !   !     do k=1,N
    !   !       if (ii(k).eq.i .and. j<jj(k)) then
    !   !         qnew2(:,i,j) = qnew(:,i,j)
    !   !       else if (ii(k).eq.i-1 .and. j<jj(k) ) then
    !   !         qnew2(:,i,j) = qnew(:,i,j)
    !   !       else if (ii(k).eq.i .and. j.gt.jj(k)+1 ) then
    !   !         qnew(:,i,j) = qold2(:,i,j)
    !   !       else if (ii(k).eq.i-1 .and. j.gt.jj(k) ) then
    !   !         qnew(:,i,j)= qold2(:,i,j)
    !   !       end if
    !   !     end do
    !   !   end do
    !   ! end do
  ! end if
  kk =1
    ! if (ids.eq.2)then
      if (kk.eq.1) then
      ! Splicing the top values and bottom values and gluing them together:
      do i=1,N
        i1 = unS_cells_i(i)
        j1 = unS_cells_j(i)
        i2 = upS_cells_i(i)
        j2 = upS_cells_j(i)
        a1 = un_area_ij(i1,j1)
        a2 = up_area_ij(i2,j2)
        if (ot) then
          qnew2(:,i1,-1:j1) = qnew(:,i1,-1:j1)
          qnew(:,i2,j2:my+2) = qnew2(:,i2,j2:my+2)
        else
          if (a1.ge.0.5d0 .and. a1.ne.1.d0 )then !(abs(a1-0.5d0).lt.1d-10 .or. a1.gt.0.5d0)) then! .and. a1.ne.1.d0) then
            qnew2(:,i1,-1:j1-1) = qnew(:,i1,-1:j1-1)  ! highlight these if overtop
          else
            qnew2(:,i1,-1:j1) = qnew(:,i1,-1:j1)
          end if
          if (a2.ge.0.5d0 .and. a2.ne.1.d0)then!(abs(a2-0.5d0).lt.1d-10 .or. a2.gt.0.5d0)) then! .and. a2.ne.1.d0) then
            qnew(:,i2,j2+1:my+2) = qnew2(:,i2,j2+1:my+2) !
          else
            qnew(:,i2,j2:my+2) = qnew2(:,i2,j2:my+2)
          end if
        end if
      end do
    ! else
    do i=1,N2
      i1 = unS_cells_i2(i)
      j1 = unS_cells_j2(i)
      i2 = upS_cells_i2(i)
      j2 = upS_cells_j2(i)
      a1 = un_area_ij(i1,j1)
      a2 = up_area_ij(i2,j2)
      if (ot) then
        qnew2(:,i1,-1:j1) = qnew(:,i1,-1:j1)
        qnew(:,i2,j2:my+2) = qnew2(:,i2,j2:my+2)
      else
        if (a1.ge.0.5d0 .and. a1.ne.1.d0 )then !(abs(a1-0.5d0).lt.1d-10 .or. a1.gt.0.5d0)) then! .and. a1.ne.1.d0) then
          qnew2(:,i1,-1:j1-1) = qnew(:,i1,-1:j1-1)  ! highlight these if overtop
        else
          qnew2(:,i1,-1:j1) = qnew(:,i1,-1:j1)
        end if
        if (a2.ge.0.5d0 .and. a2.ne.1.d0)then!(abs(a2-0.5d0).lt.1d-10 .or. a2.gt.0.5d0)) then! .and. a2.ne.1.d0) then
          qnew(:,i2,j2+1:my+2) = qnew2(:,i2,j2+1:my+2) !
        else
          qnew(:,i2,j2:my+2) = qnew2(:,i2,j2:my+2)
        end if
      end if
    end do
    end if

    ! end if
      ! qnew2=qnew

    return
    end subroutine step2ds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    subroutine rn2(ixy, maxm, num_eqn,num_waves,num_aux,num_ghost, &
             num_cells, ql, qr, auxl, auxr, fwave, s, amdq, apdq)

    ! Normal Riemann solver for the 2D SHALLOW WATER equations
    !     with topography:
    !     #        h_t + (hu)_x + (hv)_y = 0                           #
    !     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
    !     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

    ! This solver is based on David George's solver written for GeoClaw.
    ! It has been modified to be compatible with f2py (and thus PyClaw).

    ! waves:     3
    ! equations: 3

    ! Conserved quantities:
    !       1 depth
    !       2 x_momentum
    !       3 y_momentum

    ! Auxiliary fields:
    !       1 bathymetry

    ! The gravitational constant grav should be in the common block cparam.

    ! See http://www.clawpack.org/riemann.html for a detailed explanation
    ! of the Riemann solver API.

    implicit none

    real(kind=8) :: g
    real(kind=8), parameter :: drytol = 1.e-8
    !common /cparam/ grav

    integer, intent(in) :: maxm,num_eqn,num_aux,num_waves
    integer, intent(in) :: num_ghost,num_cells,ixy

    real(kind=8),intent(inout)::ql(3,maxm)
    real(kind=8),intent(inout)::qr(3,maxm)
    real(kind=8),intent(in)::auxl(maxm)
    real(kind=8),intent(in)::auxr(maxm)

    real(kind=8)::fwave(num_eqn,num_waves,maxm),s(num_waves,maxm)
    real(kind=8)::apdq(num_eqn,maxm),amdq(num_eqn,maxm)

    !local only
    integer m,i,mw,maxiter,mu,nv
    real(kind=8) wall(3)
    real(kind=8) fw(3,3)
    real(kind=8) sw(3)

    real(kind=8) hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
    real(kind=8) bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    real(kind=8) s1m,s2m
    real(kind=8) hstar,hstartest,hstarHLL,sLtest,sRtest
    real(kind=8) tw,dxdc

    logical rare1,rare2

    g = 1.d0
    m = size(ql,2)
    ! print *, "qL", ql
    ! print *, "qR", qr
    ! print* , "auxL",auxl
    ! print*, "auxR",auxr
    !loop through Riemann problems at each grid cell
    do i=2,m

    !-----------------------Initializing------------------------------
       !inform of a bad riemann problem from the start
       if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
          write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
       endif

       !Initialize Riemann problem for grid interface
       do mw=1,num_waves
            s(mw,i)=0.d0
               fwave(1,mw,i)=0.d0
               fwave(2,mw,i)=0.d0
               fwave(3,mw,i)=0.d0
       enddo

!        !set normal direction
       if (ixy.eq.1) then
          mu=2
          nv=3
       else
          mu=3
          nv=2
       endif

       !zero (small) negative values if they exist
       if (qr(1,i-1).lt.0.d0) then
             qr(1,i-1)=0.d0
             qr(2,i-1)=0.d0
             qr(3,i-1)=0.d0
       endif

       if (ql(1,i).lt.0.d0) then
             ql(1,i)=0.d0
             ql(2,i)=0.d0
             ql(3,i)=0.d0
       endif

       !skip problem if in a completely dry area
       if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
          go to 30
       endif

       !Riemann problem variables
       hL = qr(1,i-1)
       hR = ql(1,i)
       huL = qr(mu,i-1)
       huR = ql(mu,i)
       bL = auxr(i-1)
       bR = auxl(i)

       hvL=qr(nv,i-1)
       hvR=ql(nv,i)

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
          call riemanntype(hL,hL,uL,-uL,hstar,s1m,&
          s2m,rare1,rare2,1,drytol,g)
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
          call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,&
             rare1,rare2,1,drytol,g)
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

       maxiter = 1

       call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,&
          huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,&
                                      drytol,g,sw,fw)

!        !eliminate ghost fluxes for wall
       do mw=1,3
          sw(mw)=sw(mw)*wall(mw)

             fw(1,mw)=fw(1,mw)*wall(mw)
             fw(2,mw)=fw(2,mw)*wall(mw)
             fw(3,mw)=fw(3,mw)*wall(mw)
       enddo

       do mw=1,num_waves
          s(mw,i)=sw(mw)
          fwave(1,mw,i)=fw(1,mw)
          fwave(mu,mw,i)=fw(2,mw)
          fwave(nv,mw,i)=fw(3,mw)
       enddo

30      continue
    enddo



!============= compute fluctuations=============================================
       amdq(1:3,:) = 0.d0
       apdq(1:3,:) = 0.d0
       do i=2,m
          do  mw=1,num_waves
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

    return
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine rotate_state(q,q_rot,n_vec,t_vec)
    !   ! n_vec is the normal direction unit vector
    !   ! t_vec is the transverse direction unit vector, OG to n_vec
    !   ! q is the Cartesian coordinate aligned state vec
    !   ! q_rot is the rotated state vec
    !   implicit none
    !   real(8) :: q(3),q_rot(3),n_vec(2),t_vec(2)
    !   real(8) :: vel(2)
    !
    !   ! if (abs((n_vec(1)**2 + n_vec(2)**2)-1).gt.1d-8) then
    !   !   n_vec = n_vec/sqrt((n_vec(1)**2 + n_vec(2)**2))
    !   ! end if
    !   ! if (abs((t_vec(1)**2 + t_vec(2)**2)-1).gt.1d-8) then
    !   !   t_vec = t_vec/sqrt((t_vec(1)**2 + t_vec(2)**2))
    !   ! end if
    !   q_rot(1) = q(1)
    !   vel = q(2:3)
    !   q_rot(2) = vel(1)*n_vec(1) + vel(2)*n_vec(2)
    !   q_rot(3) = vel(1)*t_vec(1) + vel(2)*t_vec(2)
    ! end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine riemann_solver(ql,qr,bL,bR,s,fwave,amdq,apdq,ixy)
        implicit none
        integer :: ixy
        real(8) :: ql(3),qr(3),bL,bR,s(3),fwave(3,3),amdq(3),apdq(3)
        ! local
        real(8) :: swl(3),fwl(3,3)
        real(8)::hL,hR,uL,um,uR,huL,huR,vL,vR,hvL,hvR,wall(3),phiL,phiR
        integer :: mw,mu,nv
        real(8) :: hstar,s1m,s2m,hstartest,sL,sR,uhat,chat,sRoe1,sRoe2
        real(8) :: sE1,sE2,grav,dry_tolerance,g,drytol
        logical :: rare1,rare2
        common /cparam/ grav,dry_tolerance

        g= grav
        drytol = dry_tolerance
        print *, g,drytol
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
           call riemanntype(hL,hL,uL,-uL,hstar,um,s1m,s2m,rare1,&
                      rare2,1,drytol,g)
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
           call riemanntype(hR,hR,-uR,uR,hstar,um,s1m,s2m,rare1,&
                     rare2,1,drytol,g)
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

    end subroutine

    FUNCTION area_polygon(x, y) RESULT(fn_val)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-07-04  Time: 12:24:06

    IMPLICIT NONE

    REAL(8), INTENT(IN)     :: x(:)
    REAL(8), INTENT(IN)     :: y(:)
    INTEGER  :: nb
    REAL(8)                 :: fn_val, v1(2),v2(2)

    !*****************************************************************

    !   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
    ! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
    ! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
    ! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
    ! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
    ! itself.

    !*****************************************************************

    INTEGER  :: i, n, nm1
    REAL     :: a

    nb = size(x)
    n = nb
    a = 0.d0
    do i=1,nb-1
      v1 = (/x(i),y(i)/)
      v2 = (/x(i+1),y(i+1)/)
      a = a+v1(1)*v2(2) - v2(1)*v1(2)
    end do
    fn_val = abs(a/2.d0)
    end function

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
                                            1,drytol,g)
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

      end subroutine !r

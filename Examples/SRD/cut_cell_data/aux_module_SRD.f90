! CODES FOR SRD:


module aux_module_SRD

  use hbox_layer
  ! use redistribute2D
  implicit none
  ! real(8), parameter:: tol = 1d-7, pi=3.14159265359d0

contains
  ! function dist(A,B) result(distance)
  !   implicit none
  !   real(kind=8) :: A(2),B(2)
  !   real(kind=8) :: distance
  !   distance = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
  ! end function dist
  ! FUNCTION area_polygon(x, y) RESULT(fn_val)
  !
  ! ! Code converted using TO_F90 by Alan Miller
  ! ! Date: 2000-07-04  Time: 12:24:06
  !
  ! IMPLICIT NONE
  !
  ! REAL(8), INTENT(IN)     :: x(:)
  ! REAL(8), INTENT(IN)     :: y(:)
  ! INTEGER  :: nb
  ! REAL(8)                 :: fn_val, v1(2),v2(2)
  !
  ! !*****************************************************************
  !
  ! !   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
  ! ! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
  ! ! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
  ! ! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
  ! ! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
  ! ! itself.
  !
  ! !*****************************************************************
  !
  ! INTEGER  :: i, n, nm1
  ! REAL     :: a
  !
  ! nb = size(x)
  ! n = nb
  ! a = 0.d0
  ! do i=1,nb-1
  !   v1 = (/x(i),y(i)/)
  !   v2 = (/x(i+1),y(i+1)/)
  !   a = a+v1(1)*v2(2) - v2(1)*v1(2)
  ! end do
  ! fn_val = abs(a/2.d0)
  ! end function

  subroutine cut_cells_find(x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,mx,my,ii,jj,intersections,N_cells)   !  O
    ! returns ii, jj which are 4*mx (long) arrays that contain i-index and j-index respectively
    ! of cut cells, and N which is the actual number of cut cells (so go from ii(1:N),jj(1:N))
    implicit none

    real(8) :: x_0,y_0,x_e,y_e,dx,dy,xlower,ylower
    integer :: mx,my,ii(4*mx),jj(4*mx),N_cells,N
    real(8) :: intersections(2,4*mx)

    ! local
    integer :: i,j,k,i_0,j_0,size1,mbc
    integer, allocatable :: indexes(:)
    real(8) :: xe(-2:mx+2),ye(-2:my+2),slope_bar,theta,dist_x,dist_y
    real(8), allocatable :: intersect_0(:,:), intersect_top(:,:),midpoints(:,:)
    xe = (/(xlower+i*dx, i=-2,mx+2)/)
    ye = (/(ylower+i*dy, i=-2,my+2)/)
    mbc =2
    ! wall parameters:
    if (x_e-x_0 .ne. 0.d0) then
      slope_bar = (y_e-y_0)/(x_e-x_0)
      theta = atan(slope_bar)
    else
      i=-2
      j=-2
      do while (x_0 > xe(i))
        i = i +1
      end do
      do while  (y_0 > ye(j))
        j=j+1
      end do
      ii = i-1
      jj = (/(j-1+k,k=1,my)/)
      N_cells = my
      return
    end if

    ! get intersections between barrier and grid lines :
    do k=-2,mx+2
      if (xe(k)-tol .le. x_e .and. xe(k)+tol .ge. x_0) then
        dist_x = x_e - xe(k)
        dist_y = dist_x * tan(pi-theta)
        call AddToList_verts(intersect_0,(/xe(k),y_e+dist_y/))
      end if
    end do

    do k=-2,my+2
      if (ye(k)-tol.le.max(y_e,y_0) .and. ye(k)+tol.ge.min(y_0,y_e)) then
        dist_y = y_0 - ye(k)
        dist_x = dist_y/tan(pi-theta)
        call AddToList_verts(intersect_0,(/x_0+dist_x,ye(k)/))
      end if
    end do
    ! treat array (get rid of duplicate and sort)
       intersect_top = remove_dups(intersect_0)
       ! print*, "Intersections: x", intersect_top(1,:)
       ! print*, "intersections: y", intersect_top(2,:)
       size1 = size(intersect_top,2)
       N=size1
      allocate(indexes(size1))
      call KB07AD(intersect_top(1,:),size1,indexes)
      intersect_top(2,:) = intersect_top(2,indexes)
      ! print*,"Intersections; ", intersect_top(1,:)
      ! print*,"Intersections: ", intersect_top(2,:)
      intersections(:,1:N) = intersect_top(:,1:N)
      ! for each pair of neighboring intersections, find which cell it cuts
         !find the midpoint of the intersections and see which cell they fit in
         ! start with i_0,j_0 and walk up or down
        allocate(midpoints(N-1,2))
        do i=1, N-1
          midpoints(i,1:2) = (/0.5d0*(intersect_top(1,i)+intersect_top(1,i+1)),&
          0.5d0*(intersect_top(2,i)+intersect_top(2,i+1))/)
        end do
        ! write(*,*) "midpoint ", midpoints

      ! put them into ii and jj
      ! the first cell where the barrier starts is always cut:
      do i=1-mbc,mx+mbc
        if (midpoints(1,1) .gt. xe(i) .and. midpoints(1,1) .lt. xe(i+1)) then
          i_0 = i
          exit
        end if
      end do
      do i = 1-mbc,my+mbc
        if (midpoints(1,2) .gt. ye(i) .and. midpoints(1,2) .lt. ye(i+1)) then
          j_0 = i
          exit
        end if
      end do

      ii(1) = i_0
      jj(1) = j_0
      j = 1 ! counter for cut cells
      do i = 1, N-2
        if (midpoints(i+1,1) .lt. xe(ii(i)) .and. midpoints(i+1,1) .gt. xe(ii(i)-1)) then
          ii(i+1) = ii(i) - 1
        else if (midpoints(i+1,1) .gt. xe(ii(i)+1) .and. midpoints(i+1,1) .lt. xe(ii(i)+2)) then
          ii(i+1) = ii(i) + 1
        else if (midpoints(i+1,1) .gt. xe(ii(i)) .and. midpoints(i+1,1) .lt. xe(ii(i)+1)) then
          ii(i+1) = ii(i)
        end if
        if (midpoints(i+1,2) .lt. ye(jj(i)) .and. midpoints(i+1,2) .gt. ye(jj(i)-1)) then
          jj(i+1) = jj(i) - 1
        else if (midpoints(i+1,2) .gt. ye(jj(i)+1) .and. midpoints(i+1,2) .lt. ye(jj(i)+2)) then
          jj(i+1) = jj(i) + 1
        else if (midpoints(i+1,2) .gt. ye(jj(i)) .and. midpoints(i+1,2) .lt. ye(jj(i)+1)) then
          jj(i+1) = jj(i)
        end if
      end do
      ! N is number of intersections, N_cells number of cut cells.
      N_cells=N-1
      ! fortran indexing
      ii = ii + 1
      jj = jj + 1
   end subroutine


    subroutine small_cells_geom(ii,jj,intersections,mx,my,dx,dy,xlower,ylower,N_cells,&        ! O
      type_supper,type_sunder,lengths_supper,lengths_sunder,area_supper,area_sunder)

      implicit none
      integer :: ii(:),jj(:),mx,my,N_cells,i,i_0,j_0,j,k1,k2
      real(8) :: intersections(:,:),xlower,ylower,x1,y1,x2,y2,m,dx,dy
      real(8) :: x_0,x_e,y_0,y_e,y_int
      integer :: type_supper(4*mx), type_sunder(4*mx)
      real(8) :: area_supper(4*mx), area_sunder(4*mx)
      real(8) :: lengths_supper(5,4*mx), lengths_sunder(5,4*mx)
      real(8) :: coord_supper(2,6,4*mx), coord_sunder(2,6,4*mx)
      real(8) :: coords(2,4),  xe(-2:mx+2),ye(-2:my+2)
      integer :: uo_array(4)

      !initialization:
      lengths_supper = 0.d0
      lengths_sunder = 0.d0
      uo_array = (/3,4,1,2/) ! the order of cells' vertices checking for under-small cells' vertices

      ! edges:
      xe = (/(xlower+i*dx, i=-2,mx+2)/)
      ye = (/(ylower+i*dy, i=-2,my+2)/)

      do i=1,N_cells
        k1 = 1 ! counter for num sides for upper small cell
        k2 = 1 ! counter for num sides for under small cell
        ! the index
        i_0 = ii(i)
        j_0 = jj(i)
        ! the box
        x1 = xe(i_0-1)
        x2 = xe(i_0)
        y1 = ye(j_0-1)
        y2 = ye(j_0)
        coords(1,:) = (/x1,x1,x2,x2/)
        coords(2,:) = (/y1,y2,y2,y1/)
        ! the bar
        x_0 = intersections(1,i)
        y_0 = intersections(2,i)
        x_e = intersections(1,i+1)
        y_e = intersections(2,i+1)
        m = (y_e-y_0)/(x_e-x_0)
        y_int = y_0-m*x_0
        ! the coordinates of small cells (upper and under)
        coord_supper(:,1,i) = (/x_0,y_0/)
        coord_supper(:,2,i) = (/x_e,y_e/)
        k1 = k1+2
        coord_sunder(:,1,i) = (/x_0,y_0/)
        coord_sunder(:,2,i) = (/x_e,y_e/)
        k2 = k2+2
        do j=1,4
          if (m*coords(1,4-(j-1))+y_int < coords(2,4-(j-1))) then
            coord_supper(:,k1,i) = coords(:,4-(j-1))
            k1= k1 + 1
          end if
          if (m*coords(1,uo_array(j))+y_int > coords(2,uo_array(j))) then
            coord_sunder(:,k2,i) = coords(:,uo_array(j))
            k2 = k2 + 1
          end if
        end do
        coord_supper(:,k1,i) = (/x_0,y_0/)
        coord_sunder(:,k2,i) = (/x_0,y_0/)

        ! side lengths
        do j=1,k1-1
          lengths_supper(j,i) = dist(coord_supper(:,j,i),coord_supper(:,j+1,i))/dx
        end do
        do j=1,k2-1
          lengths_sunder(j,i) = dist(coord_sunder(:,j,i),coord_sunder(:,j+1,i))/dx
        end do
        k1 = k1-1 ! actual number of sides
        k2 = k2-1
        ! type 1-8
        if (k1.eq.5 .and. m>0) then
          type_supper(i) = 1
          type_sunder(i) = 1
        else if (k1.eq.4 .and. m>0) then
          if (abs(m)<1) then
            type_supper(i) = 3
            type_sunder(i) = 3
          else
            type_supper(i) = 2
            type_sunder(i) = 2
          endif
        else if (k1.eq.3 .and. m>0) then
          type_supper(i) = 4
          type_sunder(i) = 4
        else if (k1.eq.5 .and. m<0) then
          type_supper(i) = 8
          type_sunder(i) = 8
        else if (k1.eq.4 .and. m<0) then
          if (abs(m)<1) then
            type_supper(i) = 6
            type_sunder(i) = 6
          else
            type_supper(i) = 7
            type_sunder(i) = 7
          end if
        else if (k1.eq.3 .and. m<0) then
          type_supper(i) = 5
          type_sunder(i) = 5
        end if
        ! area of small cells
        area_supper(i) = area_polygon(coord_supper(1,1:k1+1,i), coord_supper(2,1:k1+1,i))/(dx*dy)
        area_sunder(i) = area_polygon(coord_sunder(1,1:k2+1,i), coord_sunder(2,1:k2+1,i))/(dx*dy)
      end do


    end subroutine


  !   subroutine small_cells_update_correct(ii,jj,dtdx,qold,qold2,qnew,qnew2,aux,&              !  [   ]
  !     intersections,type_supper,type_sunder,lengths_supper,lengths_sunder,&
  !       area_sunder,area_supper,ixy,N_cells,xlower,ylower,xupper,yupper,&
  !       mx,my,mbc,x_0,y_0,x_e,y_e,wall_height,fm,fp,gm,gp,fm2,fp2,gm2,gp2)
  !     ! takes qnew/qnew2 and uses qold/qold2 to correct the update so that its specific to the small cells' geometry
  !     ! ie computes Q^ in the paper for both up and under small cells
  !     implicit none
  !
  !     ! input / output
  !     integer :: ii(:),jj(:),type_supper(:),type_sunder(:),ixy,N_cells,mx,my,mbc
  !     real(8), intent(in) :: qold(3,1-mbc:mx+mbc,1-mbc:my+mbc),qold2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) :: qnew(3,1-mbc:mx+mbc,1-mbc:my+mbc),qnew2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) :: aux(1,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) :: intersections(:,:),lengths_sunder(:,:),lengths_supper(:,:)
  !     real(8) :: area_supper(:),area_sunder(:),xlower,ylower,xupper,yupper
  !     real(8) :: x_0,y_0,x_e,y_e,wall_height,dtdx
  !     real(8) fp(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) fm(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) gp(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) gm(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) fp2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) fm2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) gp2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !     real(8) gm2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
  !
  !     ! local
  !     real(8) :: n_vec(2),t_vec(2),x,y
  !     real(8) :: q_hbox_d(3,mx*4), q_hbox_u(3,mx*4), aux_hbox_u(mx*4),aux_hbox_d(mx*4)
  !     logical :: L2R, R2L
  !     integer :: is,js,i,j,k,m,num_frags_d(mx*4),index_frags_d(2,4,mx*4)
  !     integer :: num_frags_u(mx*4),index_frags_u(2,4,mx*4),look(2)
  !     real(8) :: hbox_areas_u(mx*4), area_frags_u(4,mx*4)
  !     real(8) :: hbox_areas_d(mx*4), area_frags_d(4,mx*4)
  !     real(8) :: wave_wall(3,3), amdq_wall(3),apdq_wall(3),s_wall(3)
  !     real(8) :: hL,hR,huL,huR,hvL,hvR,bL,bR,hstarL,hstarR,ustarL,ustarR
  !     logical :: lexist
  !
  !     ! get the hboxes for normal barrier flux calculation:
  !     inquire (file="./hbox_data.txt",exist=lexist)
  !     ! print *, "EXIST:",lexist
  !     open (unit=2,file="./hbox_data.txt",FORM="FORMATTED",STATUS="OLD",&
  !     ACTION="READ",access='sequential')
  !     rewind 2
  !     read(2,*) m
  !     read(2,*) num_frags_d
  !     read(2,*) num_frags_u
  !     read(2,*) index_frags_d
  !     read(2,*) index_frags_u
  !     read(2,*) hbox_areas_d
  !     read(2,*) hbox_areas_u
  !     read(2,*) area_frags_u
  !     read(2,*,end=100) area_frags_d
  !     close(2)!,status="keep")
  ! 100 continue
  !     do i = 1,m
  !       do j=1,num_frags_d(i)
  !         look = index_frags_d(:,j,i)
  !         q_hbox_d(:,i) = q_hbox_d(:,i) + qold(:,look(1),look(2))*area_frags_d(j,i)
  !         aux_hbox_d(i) = aux_hbox_d(i) + aux(1,look(1),look(2))*area_frags_d(j,i)
  !       end do
  !     enddo
  !     do i = 1,m
  !       do j=1,num_frags_u(i)
  !         look = index_frags_u(:,j,i)
  !         q_hbox_u(:,i) = q_hbox_u(:,i) + qold2(:,look(1),look(2))*area_frags_u(j,i)
  !         aux_hbox_u(i) = aux_hbox_u(i) + aux(1,look(1),look(2))*area_frags_u(j,i)
  !       end do
  !     enddo
  !     ! call down_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii(1),jj(1),x_0,y_0,&
  !     ! ii(N_cells),jj(N_cells),x_e,y_e,qold,aux,1,q_hbox_d,aux_hbox_d)
  !     ! call up_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii(1),jj(1),x_0,y_0,&
  !     ! ii(N_cells),jj(N_cells),x_e,y_e,qold2,aux,1,q_hbox_u,aux_hbox_u)
  !     ! for x swipe
  !     print *, "HBOX HTS DOWN: ------------------------"
  !     print*,   q_hbox_d(1,1:m)
  !     print *, "---------------------------------------"
  !     print *, "HBOX HTS UP: --------------------------"
  !     print *, q_hbox_u(1,1:m)
  !     print *, "----------------------------------------"
  !     if (ixy.eq.1) then
  !
  !       do i = 1,N_cells
  !         is = ii(i)
  !         js = jj(i)
  !         x = intersections(1,i+1)- intersections(1,i)
  !         y = intersections(2,i+1) - intersections(2,i)
  !         n_vec = (/-y,x/)
  !         n_vec = n_vec/(sqrt(x**2+y**2))
  !         t_vec = (/x,y/)
  !         t_vec = t_vec/(sqrt(x**2+y**2))
  !
  !         call rotate_state(q_hbox_d(:,i),q_hbox_d(:,i), &
  !         n_vec,t_vec)
  !         call rotate_state(q_hbox_u(:,i),q_hbox_u(:,i),&
  !          n_vec,t_vec)
  !
  !         hL = q_hbox_u(1,i)
  !         hR = q_hbox_d(1,i)
  !         huL = q_hbox_u(2,i)
  !         huR = q_hbox_d(2,i)
  !         hvL= q_hbox_u(3,i)
  !         hvR = q_hbox_d(3,i)
  !         bL = aux_hbox_u(i)
  !         bR = aux_hbox_d(i)
  !
  !         call barrier_passing(hL,hR,huL,huR,bL,bR,wall_height,&
  !                   L2R,R2L,hstarL,hstarR,ustarL,ustarR)
  !         call redistribute_fwave(1,q_hbox_u(:,i),q_hbox_d(:,i),aux_hbox_u(i),&
  !            aux_hbox_d(i),wall_height,1,wave_wall,s_wall,amdq_wall,apdq_wall,3,&
  !            3,L2R,R2L)
  !         ! print *, hL, hR
  !         ! print *,"amdq", amdq_wall
  !         ! print*, "apdq", apdq_wall
  !
  !
  !         ! print*, "*******TYPE: ", type_supper(i), "***********"
  !         ! print*, "the undersmall cell:" , qnew(:,is,js)
  !         ! print*, "the uppersmall cell", qnew2(:,is,js)
  !
  !         select case (type_supper(i))
  !         case (1) ! TYPE 1 Cut cell
  !           ! upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*lengths_supper(1,i)*&
  !             (amdq_wall) !f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*(lengths_supper(2,i))*fm2(:,is+1,js)!&
  !           !  + f(qold2(:,is,js),1))
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*lengths_sunder(1,i)*&
  !             ( apdq_wall )!f(qold(:,is,js),1)
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*lengths_sunder(2,i)*(fm(:,is+1,js))
  !           !   +f(qold(:,is,js),1)&
  !         case (2) ! TYPE 2 Cut cell
  !           ! upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*lengths_supper(1,i)*&
  !             ( amdq_wall)!f(qold2(:,is,js),1)+
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*fp2(:,is,js)
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*lengths_sunder(1,i)*&
  !             (apdq_wall) !f(qold(:,is,js),1)
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*fm(:,is+1,js)
  !         case (3,6) ! TYPE 3/6 Cut cell
  !           ! upper small cell
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*((lengths_supper(2,i))*fm2(:,is+1,js) + &
  !             (lengths_supper(4,i))*fp2(:,is,js)) ! + dtdx* (lengths_supper(2,i)-lengths_supper(4,i))&
  !               ! * (f(qold2(:,is,js),1))
  !           ! under small cell
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*((lengths_sunder(2,i))*fm(:,is+1,js) + &
  !             (lengths_sunder(4,i))*fp(:,is,js))!+ dtdx* (lengths_sunder(2,i)-lengths_sunder(4,i))&
  !               ! * (f(qold(:,is,js),1))
  !
  !         case (4) ! TYPE 4 Cut cell
  !           ! upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*lengths_supper(1,i)*&
  !             ( amdq_wall)!f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*lengths_supper(3,i)*(fp2(:,is,js))
  !           ! f(qold2(:,is,js),1)-&
  !
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*lengths_sunder(1,i)*&
  !             (apdq_wall)!f(qold(:,is,js),1) +
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*(lengths_sunder(5,i))*(fp(:,is,js))
  !           ! -f(qold(:,is,js),1)&
  !
  !         case (5) ! TYPE 5 Cut cell
  !           ! upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*&!lengths_supper(1,i)*&
  !             (amdq_wall)!f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*lengths_supper(2,i)*(fm2(:,is+1,js))
  !           ! f(qold2(:,is,js),1) +
  !           ! under small cell
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*&!lengths_sunder(1,i)*&
  !             (-1*apdq_wall)!f(qold(:,is,js),1) -
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*(1-lengths_sunder(2,i))*(fm(:,is+1,js))
  !           ! f(qold(:,is,  js),1) +
  !         case (7)
  !           !  upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*&!lengths_supper(1,i)*&
  !             (amdq_wall) !f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) + dtdx*fm2(:,is+1,js)
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*&!lengths_sunder(1,i)*&
  !             (-1* apdq_wall)!f(qold(:,is,js),1)
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) + dtdx*fp(:,is,js)
  !         case(8)
  !           ! upper small cell
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*&!lengths_supper(1,i)*&
  !             (amdq_wall)!f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*(1-lengths_supper(5,i))*(-1*fp2(:,is,js))
  !           ! f(qold2(:,is,&  js),1)
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*&!lengths_sunder(1,i)*&
  !             (-1* apdq_wall)!f(qold(:,is,js),1)
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) + dtdx*lengths_sunder(3,i)*(-1*fp(:,is,js))
  !           ! f(qold(:,is,js),1)&
  !
  !         end select
  !         ! print*, "the undersmall cell: post" , qnew(:,is,js)
  !         ! print*, "the uppersmall cell: post", qnew2(:,is,js)
  !       end do
  !     end if
  !
  !    ! for y swipe
  !     if (ixy.eq.2) then
  !
  !       do i = 1,N_cells
  !         is = ii(i)
  !         js = jj(i)
  !
  !         x = intersections(1,i+1)- intersections(1,i)
  !         y = intersections(2,i+1) - intersections(2,i)
  !         n_vec = (/-y,x/)
  !         n_vec = n_vec/(sqrt(x**2+y**2))
  !         t_vec = (/x,y/)
  !         t_vec = t_vec/(sqrt(x**2+y**2))
  !
  !         call rotate_state(q_hbox_d(:,i),q_hbox_d(:,i), &
  !         n_vec,t_vec)
  !         call rotate_state(q_hbox_u(:,i),q_hbox_u(:,i),&
  !          n_vec,t_vec)
  !
  !         hL = q_hbox_u(1,i)
  !         hR = q_hbox_d(1,i)
  !         huL = q_hbox_u(2,i)
  !         huR = q_hbox_d(2,i)
  !         hvL= q_hbox_u(3,i)
  !         hvR = q_hbox_d(3,i)
  !         bL = aux_hbox_u(i)
  !         bR = aux_hbox_d(i)
  !
  !         call barrier_passing(hL,hR,huL,huR,bL,bR,wall_height,&
  !                   L2R,R2L,hstarL,hstarR,ustarL,ustarR)
  !         call redistribute_fwave(1,q_hbox_u(:,i),q_hbox_d(:,i),aux_hbox_u(i),&
  !            aux_hbox_d(i),wall_height,1,wave_wall,s_wall,amdq_wall,apdq_wall,3,&
  !            3,L2R,R2L)
  !            ! print *, hL, hR
  !            !
  !            ! print *,"amdq", amdq_wall
  !            ! print*, "apdq", apdq_wall
  !
  !        ! print*, "*******TYPE: ", type_supper(i), "***********"
  !        !
  !        !  print*, "the undersmall cell:" , qnew(:,is,js)
  !        !   print*, "the uppersmall cell", qnew2(:,is,js)
  !
  !         select case (type_supper(i))
  !         case (1) ! TYPE 1 Cut cell
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*(lengths_supper(5,i))*(gp2(:,is,js))
  !           ! f(qold2(:,is,js),2)-&
  !           !
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*lengths_sunder(3,i)*(gp(:,is,js))
  !         ! f(qold(:,is,js),2)&  -
  !
  !         case (2) ! TYPE 2 Cut cell
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js) -dtdx*((lengths_supper(2,i))*gm2(:,is,js+1)+&
  !              (lengths_supper(4,i)) * gp2(:,is,js)) !+dtdx*(lengths_supper(2,i)-lengths_supper(4,i))*&
  !                ! f(qold2(:,is,js),2)
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js) -dtdx*((lengths_sunder(2,i))*gm(:,is,js+1)+&
  !              (lengths_sunder(4,i)) * gp(:,is,js))! +dtdx*(lengths_sunder(2,i)-lengths_sunder(4,i))*&
  !                ! f(qold(:,is,js),2)
  !         case (3,6)
  !           ! upper small cell
  !           call rotate_state(qold2(:,is,js),qold2(:,is,js),n_vec,t_vec)
  !           qnew2(:,is,js) = qold2(:,is,js) - dtdx*lengths_supper(1,i)*&
  !             (amdq_wall) !f(qold2(:,is,js),1) +
  !           call rotate_state(qnew2(:,is,js),qnew2(:,is,js),t_vec,-n_vec)
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*gm2(:,is,js+1)
  !           ! under small cell
  !           call rotate_state(qold(:,is,js),qold(:,is,js),n_vec,t_vec)
  !           qnew(:,is,js) = qold(:,is,js) - dtdx*lengths_sunder(1,i)*&
  !             (apdq_wall)!f(qold(:,is,js),1)
  !           call rotate_state(qnew(:,is,js),qnew(:,is,js),t_vec,-n_vec)
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*gp(:,is,js)
  !         case (4)
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*lengths_supper(2,i)* &
  !             (gm2(:,is,js+1))!f(qold2(:,is,js),2) +
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js) - dtdx*(lengths_sunder(2,i))*&
  !             (gm(:,is,js+1))!f(qold(:,is,js),2) +
  !         case (5)
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js) - dtdx*lengths_supper(3,i)* &
  !             (gm2(:,is,js+1))!f(qold2(:,is,js),2) +
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js) -dtdx*(1-lengths_sunder(5,i))*&
  !             (gm(:,is,js+1)) !f(qold(:,is,js),2) +
  !         case (7)
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js)- dtdx* ((1-lengths_supper(4,i))*gm2(:,is,js+1) &
  !             + (1-lengths_supper(2,i))*gp2(:,is,js))!+dtdx*(lengths_supper(4,i)-lengths_supper(2,i))&
  !               ! * (f(qold2(:,is,js),2))
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js)- dtdx* ((1-lengths_sunder(4,i))*gm(:,is,js+1) &
  !             + (1-lengths_sunder(2,i))*gp(:,is,js)) !+dtdx*(lengths_sunder(4,i)-lengths_sunder(2,i))&
  !               ! * (f(qold(:,is,js),2))
  !         case(8)
  !           ! upper small cell
  !           qnew2(:,is,js) = qnew2(:,is,js) + dtdx*(1-lengths_supper(2,i))*&
  !             (-1*gp2(:,is,js)) !f(qold2(:,is,js),2)
  !           ! under small cell
  !           qnew(:,is,js) = qnew(:,is,js) + dtdx*(lengths_sunder(2,i))*&
  !             (-1*gp(:,is,js)) !f(qold(:,is,js),2)
  !         end select
  !         ! print*, "the undersmall cell: post" , qnew(:,is,js)
  !         ! print*, "the uppersmall cell: post", qnew2(:,is,js)
  !       end do
  !     end if
  !   end subroutine


    subroutine SRD_undercells(N_cells,ii,jj,mx,my,area_sunder,x_0,y_0,x_e,y_e,dx,dy,&           !  [      ]
      unS_cells_i,unS_cells_j,N_ij,all_undercells_i,all_undercells_j,k_count)
      implicit none

      integer :: mx,my,ii(:),jj(:),N_cells,N_ij(-1:mx+2,-1:my+2)
      real(8) :: area_sunder(:),x_0,y_0,x_e,y_e,m,dx,dy
      integer :: unS_cells_i(mx*4), unS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
      integer :: all_undercells_i(mx*4), all_undercells_j(mx*4) ! THESE will just give you in order the indices of all affected cells by nhood inclusions

      integer :: i,j,k_count
      logical :: TF
      ! N_ij is the matrix of number of neighborhodds each cell belongs to
      N_ij = 1 ! i.e. everybody is its own neighbor

      ! slope of barrier :
      m = (y_e-y_0)/(x_e-x_0)
      ! initialization of indices that will have to be looped over for SRD updates
      unS_cells_i = huge(1)
      unS_cells_j = huge(1)

      if (abs(m).le.1.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          else if (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i) - 1
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.lt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          elseif (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i) - 1
            unS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.gt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          elseif (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i) + 1
            unS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      end if

      k_count = 1
      do i = 1,N_cells
        all_undercells_i(k_count) = ii(i)
        all_undercells_j(k_count) = jj(i)
        k_count = k_count + 1
        if (area_sunder(i) .lt. 0.5d0 .and. abs(area_sunder(i)-0.5d0).gt.5d-13) then
          TF = check_intpair_in_set((/unS_cells_i(i),unS_cells_j(i)/),&
          ii(max(1,i-1):min(i+1,N_cells)),jj(max(1,i-1):min(i+1,N_cells)))
          if (.not. TF) then
            all_undercells_i(k_count) = unS_cells_i(i)
            all_undercells_j(k_count) = unS_cells_j(i)
            k_count = k_count + 1
          end if
        end if
      end do

      k_count = k_count - 1 ! number of actually affected cells for which SRD update applies

    end subroutine

    subroutine SRD_uppercells(N_cells,ii,jj,mx,my,area_supper,x_0,y_0,x_e,y_e,dx,dy,&        !  [     ]
      upS_cells_i,upS_cells_j,N_ij,all_uppercells_i,all_uppercells_j,k_count)
      implicit none

      integer :: mx,my,ii(:),jj(:),N_cells,N_ij(-1:mx+2,-1:my+2)
      real(8) :: area_supper(:),x_0,y_0,x_e,y_e,m,dx,dy
      integer :: upS_cells_i(mx*4), upS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
      integer :: all_uppercells_i(mx*4), all_uppercells_j(mx*4) ! THESE will just give you in order the indices of all affected cells by nhood inclusions  (AKA the YELLOW ARRAY OF INDICES)

      integer :: i,j,k_count
      logical :: TF
      ! N_ij is the matrix of number of neighborhodds each cell belongs to
      N_ij = 1

      ! slope of barrier :
      m = (y_e-y_0)/(x_e-x_0)
      ! initialization of indices that will have to be looped over for SRD updates
      upS_cells_i = huge(1)
      upS_cells_j = huge(1)

      if (abs(m).le.1.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          else if (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i) + 1
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.lt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          elseif (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i) + 1
            upS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.gt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          elseif (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i) - 1
            upS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      end if

      k_count = 1
      do i = 1,N_cells
        all_uppercells_i(k_count) = ii(i)
        all_uppercells_j(k_count) = jj(i)
        k_count = k_count + 1
        if (area_supper(i) .lt. 0.5d0.and. abs(area_supper(i)-0.5d0).gt.5d-13) then
          TF = check_intpair_in_set((/upS_cells_i(i),upS_cells_j(i)/),&
          ii(max(1,i-1):min(i+1,N_cells)),jj(max(1,i-1):min(i+1,N_cells)))
          if (.not. TF) then
            all_uppercells_i(k_count) = upS_cells_i(i)
            all_uppercells_j(k_count) = upS_cells_j(i)
            k_count = k_count + 1
          end if
        end if
      end do

      k_count = k_count - 1 ! number of actually affected cells for which SRD update applies

    end subroutine


    function check_intpair_in_set(intpair,setx,sety) result (TF)
      implicit none
      integer :: intpair(2), setx(:),sety(:),i,n
      logical :: TF
      TF = .false.
      n = size(setx)
      do i =1,n
        if (intpair(1) .eq. setx(i) .and. intpair(2) .eq. sety(i)) then
          TF = .true.
          return
        end if
      end do
    end function


    subroutine area_cells(ii,jj,mx,my,mbc,area_sunder,area_supper,up_area_ij,un_area_ij,keep)
      implicit none
      integer :: ii(:),jj(:),mx,my,mbc
      real(8) :: area_sunder(:),area_supper(:)
      real(8),intent(inout) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
      logical :: keep
      ! local
      integer :: N,i

      ! initialization of area matrix
      if (.not. keep) then
        up_area_ij = 1.d0
        un_area_ij = 1.d0
      endif

      ! change the small cell ones
      N = size(ii)
      do i=1,N
        up_area_ij(ii(i),jj(i)) = area_supper(i)
        un_area_ij(ii(i),jj(i)) = area_sunder(i)
      end do
    end subroutine


    subroutine SRD_update_correct(qnew,qnew2,all_undercells_i,all_undercells_j,&      !  [     ]
                    all_uppercells_i,all_uppercells_j,ii,jj,N_ij_up,N_ij_un,&
                    unS_cells_i,unS_cells_j,upS_cells_i,upS_cells_j,up_area_ij,&
                    un_area_ij,mx,my,mbc)
      implicit none
      real(8) :: qnew(3,1-mbc:mx+mbc,1-mbc:my+mbc),qnew2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: all_undercells_i(:),all_uppercells_i(:),all_undercells_j(:)
      integer :: all_uppercells_j(:),ii(:),jj(:),N_ij_up(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: N_ij_un(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: unS_cells_i(:),unS_cells_j(:),upS_cells_i(:),upS_cells_j(:)
      real(8) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc
      ! local
      integer :: i,N,i0,j0,num_ij,i_n,j_n, array(mx*4),k,array2(mx*4)
      real(8) :: beta,alpha,nhood_val(3)

      ! DO UNDER SIDE FIRST:
      N = size(all_undercells_i)
      k = 1
      do i=1,N
        i0 = all_undercells_i(i)
        j0 = all_undercells_j(i)
        if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
          array(i) = k
          k = k +1
        else
          array(i)=k-1
        end if
        ! print *, " I0,J0:", i0,j0
        ! print *, "ARRAY:", array(i)
        num_ij = N_ij_un(i0,j0)
        ! print *, "NUM_IJ:", num_ij
        if (num_ij .eq. 2) then
          qnew(:,i0,j0) = qnew(:,i0,j0)/2.d0
        else if (num_ij .eq. 1) then
          ! get neighborhood value:
          ! print *, "area: UN", un_area_ij(i0,j0)
          if (un_area_ij(i0,j0).lt.0.5d0 .and. abs(un_area_ij(i0,j0)-0.5d0).gt.5d-13)then
            i_n = unS_cells_i(array(i))
            j_n = unS_cells_j(array(i))
            ! print*, "IN,JN:", i_n,j_n
            beta = un_area_ij(i_n,j_n)
            alpha = un_area_ij(i0,j0)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew(:,i0,j0)+&
                 beta/2.d0*qnew(:,i_n,j_n))
            qnew(:,i0,j0) = nhood_val
            qnew(:,i_n,j_n) = qnew(:,i_n,j_n) + nhood_val/2.d0
          end if
        end if
      end do

      ! DO UPPER SIDE FIRST:
      N = size(all_uppercells_i)
      k=1
      do i=1,N
        i0 = all_uppercells_i(i)
        j0 = all_uppercells_j(i)
        if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
          array2(i) = k
          k =k +1
        else
          array2(i)=k-1
        end if
        ! print*, "ARRAY2:", array2(i)
        num_ij = N_ij_up(i0,j0)
        ! print *, "NUM_IJ:", num_ij
        if (num_ij .eq. 2) then
          qnew2(:,i0,j0) = qnew2(:,i0,j0)/2.d0
        else if (num_ij .eq. 1) then
          ! print *, "area: UP", up_area_ij(i0,j0)
          ! get neighborhood value:
          if (up_area_ij(i0,j0).lt.0.5d0 .and. abs(up_area_ij(i0,j0)-0.5d0).gt.5d-13)then
            i_n = upS_cells_i(array2(i))
            j_n = upS_cells_j(array2(i))
            ! print *, "In,jn:",i_n,j_n
            beta = up_area_ij(i_n,j_n)
            alpha = up_area_ij(i0,j0)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew2(:,i0,j0)+&
                 beta/2.d0*qnew2(:,i_n,j_n))
            qnew2(:,i0,j0) = nhood_val
            qnew2(:,i_n,j_n) = qnew2(:,i_n,j_n) + nhood_val/2.d0
          end if
        end if
      end do


    end subroutine

    subroutine rotate_state(q,q_rot,n_vec,t_vec)
      ! n_vec is the normal direction unit vector
      ! t_vec is the transverse direction unit vector, OG to n_vec
      ! q is the Cartesian coordinate aligned state vec
      ! q_rot is the rotated state vec
      implicit none
      real(8) :: q(3),q_rot(3),n_vec(2),t_vec(2)
      real(8) :: vel(2)

      ! if (abs((n_vec(1)**2 + n_vec(2)**2)-1).gt.1d-8) then
      !   n_vec = n_vec/sqrt((n_vec(1)**2 + n_vec(2)**2))
      ! end if
      ! if (abs((t_vec(1)**2 + t_vec(2)**2)-1).gt.1d-8) then
      !   t_vec = t_vec/sqrt((t_vec(1)**2 + t_vec(2)**2))
      ! end if
      q_rot(1) = q(1)
      vel = q(2:3)
      q_rot(2) = vel(1)*n_vec(1) + vel(2)*n_vec(2)
      q_rot(3) = vel(1)*t_vec(1) + vel(2)*t_vec(2)
    end subroutine










end module

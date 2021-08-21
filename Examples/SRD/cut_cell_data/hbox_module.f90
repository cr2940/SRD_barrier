module hbox_layer
  use aux_module_hbox
  implicit none
  real(8), parameter:: tol = 5d-13, pi=3.14159265359d0

contains

  subroutine down_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,i_0,j_0,x_0,y_0,i_e,j_e,x_e,y_e,&
    m,hbox_areas,num_frags,area_frags,index_frags,layer)!,ii,jj)
    ! make it output first how many hboxes (m :: int), then for each hbox, how many fragments (num_frags), what is the area of each fragment (area_frags),
    ! what is the index loc of each fragment (index_frags), waht is the hbox total area (hbox_areas)
    implicit none
    ! input:
    real(8):: xlower,ylower,xupper,yupper,x_0,y_0,x_e,y_e ! xlower,xupper,ylower,yupper is the grid boundary points; x_0,y_0,x_e,y_e endpoints of barrier segment
    integer :: mbc,mx,my,i_0,j_0,i_e,j_e,layer !mbc,mx,my are num ghost, grid cells; i_0,j_0,i_e,j_e indices of endpoints of barrier
                                               !layer is integer indicating first or second layer of hboxes
    real(8) :: q_old(3,1-mbc:mx+mbc,1-mbc:my+mbc) ! q_old is q values on the grid at time t^n, but this should be q_old_lower as opposed to q_old_upper
                                                  ! q_old_lower is the grid q values where the cut cells are replaced with the bottom small cell values
                                                  ! q_old_upper is the grid q values where the cut cells are replaced with the top small cell values
                                                  ! this is to simplify getting the hbox averages by indexing in the context of grid
    real(8) :: aux(1,1-mbc:mx+mbc,1-mbc:my+mbc)
    ! output:
    integer :: m  ! m is number of hboxes
    integer :: num_frags(mx*4), index_frags(2,4,mx*4) ! number of fragments in each hbox ; the index location of each fragment (<=4)
    real(8) :: hbox_areas(mx*4), area_frags(4,mx*4)  ! area of each hbox ; areas (<=4) of each fragment in each hbox
    ! integer,allocatable :: ii(:),jj(:)
    ! local:
    real(8) :: dx, dy, xe(-mbc:mx+mbc), ye(-mbc:my+mbc),x1,x2,y1,y2
    real(8) :: slope_bar,theta,dist_x,dist_y,dist_x2,dist_y2,dist_x3
    real(8) :: x_1,y_1, xrange(2),yrange(2),angle,xrange2(2),yrange2(2)
    real(8) :: cands(2,4),cen_x,cen_y,A(2),B(2),C(2),D(2),sides(4,2,2),P(2)
    real(8) :: x_vers(4), y_vers(4),xmin,xmax,ymin,ymax,normal_slope, c_int
    real(8) :: last_ver(2), pen_ver(2),x_phys_ver,y_phys_ver,phiz1(2),phiz(2)
    real(8) :: d_check,x_phiz,y_phiz,v1(2),v2(2),normal_slope1,normal_slope2
    real(8) :: normal_slope3,normal_slope4,c1,c2,c3,c4,phiz2(2),x_phiz1,x_phiz2
    real(8) :: y_phiz1,y_phiz2,x_cat(2,4),y_cat(2,4),ver_x(2),v3(2),v4(2),v5(2),v6(2)
    real(8) :: e1(2,2),hbox_area,init_area,centroid(2),h_edge(2,2),last_elt(2,2),z
    real(8) :: v_last1(2),v_last2(2),cand(2),ver3(2),ver2(2),ver1(2),vec1(2),vec2(2)
    real(8) :: area,v_edge(2,2),dist_y3
    real(8),allocatable ::intersect_0(:,:), pair_inter_norm(:,:),pair_inter_par1(:,:)
    real(8),allocatable :: intersect_n(:,:),intersect_p(:,:),pair_inter_par2(:,:)
    real(8),allocatable :: intersect_top(:,:), ver_hbox(:,:),ver_hbox2(:,:),phizverts(:,:)
    real(8),allocatable :: we(:,:), wc(:,:), we_1(:,:), wc_1(:,:),hbox_cen(:,:),x_list(:)
    real(8),allocatable :: pair_inter_n(:,:), pair_inter_p1(:,:),pair_inter_p2(:,:)
    real(8),allocatable :: edges(:,:,:),inter_chk_n(:,:),inter_chk(:,:),inter_chkd(:,:)
    real(8),allocatable :: edges_chk(:,:,:),areas(:), cyclea(:,:),horz_edges(:,:,:)
    real(8),allocatable :: vert_edges(:,:,:),region_d(:,:,:),phys_ver(:,:)

    integer :: i,i2,j,j2,k1,k2,k3,i_0_p,i_0_t,j_0_p,j_0_t,k,ci,i0,r0,size1,num_edges
    integer :: num_norm, num_par,p1,m1,m2,m3
    integer:: i_1,j_1,in1,in2,i1,bot_guess(2),bot_guess2(2),i_h,j_h,chk_in,code
    integer :: counter,counter2,num_edge,ord(2), right_ind(2),isign,i3,look(2)
    integer :: where_sum(4),ml,nl,vl,wl,dummy(4)
    integer,allocatable :: rough_ind(:,:),wherea(:,:),listt(:),find(:),o(:)
    integer , allocatable :: indexes(:), odd(:),even(:),i_chk(:),j_chk(:),i_range(:),j_range(:)
    integer, allocatable :: no_int_norm(:), no_int_par(:),i_chk2(:),j_chk2(:),i_range2(:),j_range2(:)
    integer , allocatable :: no_int_norm2(:),no_edges(:),inds(:,:)
    integer , allocatable :: phys_ver_in(:),min_loc(:),max_loc(:),ver_min(:),ver_max(:)

    logical :: vertical,yes_or_no

    area_frags=0.d0
    hbox_areas = 0.d0

    dx = (xupper-xlower)/mx
    dy = (yupper-ylower)/my
    xe = (/(xlower+i*dx, i=-2,mx+2)/)
    ye = (/(ylower+i*dy,i=-2,my+2)/)

    i_0_p = i_0
    i_0_t = i_0
    j_0_p = j_0
    j_0_t = j_0

    ! wall parameters:
    vertical = .false.
    if (x_e-x_0 .ne. 0.d0) then
      slope_bar = (y_e-y_0)/(x_e-x_0)
      theta = atan(slope_bar)
    else
      vertical = .true.
    end if

    ! get intersections between barrier and grid lines :
    do k=-2,mx+2
      if (xe(k) .le. x_e+tol .and. xe(k) .ge. x_0-tol) then
        dist_x = x_e - xe(k)
        dist_y = dist_x * tan(pi-theta)
        call AddToList_verts(intersect_0,(/xe(k),y_e+dist_y/))
      end if
    end do

    do k=-2,my+2
      if (ye(k).le.max(y_e,y_0)+tol .and. ye(k).ge.min(y_0,y_e)-tol) then
        dist_y = y_0 - ye(k)
        dist_x = dist_y/tan(pi-theta)
        call AddToList_verts(intersect_0,(/x_0+dist_x,ye(k)/))
      end if
    end do

    ! treat array (get rd of duplicate and sort)
       intersect_top = remove_dups(intersect_0)
       ! print*, "Intersections: x", intersect_top(1,:)
       ! print*, "intersections: y", intersect_top(2,:)
       size1 = size(intersect_top,2)
      allocate(indexes(size1))
      call KB07AD(intersect_top(1,:),size1,indexes)
      intersect_top(2,:) = intersect_top(2,indexes)
      ! print*, "Intersections: x", intersect_top(1,:)
      ! print*, "intersections: y", intersect_top(2,:)

      ! first layer of hbox vertices
      allocate(ver_hbox(2,2*size1))
      allocate(ver_hbox2(2,2*size1))

      allocate(odd(size1))
      allocate(even(size1))
      odd = (/ (i,i=1,2*size1,2)/)
      even = (/(i,i=2,2*size1,2)/)
      ver_hbox(:,odd) = intersect_top
      do i = 1,size1
        ver_hbox(:,even(i)) = intersect_top(:,i) + (/dx*cos(pi/2.d0-theta),-dy*sin(pi/2.d0-theta)/)
      end do

      ver_hbox2(:,odd) = ver_hbox(:,even)
      do i= 1,size1
        ver_hbox2(:,even(i)) = intersect_top(:,i) + (/2*dx*cos(pi/2.d0-theta),-2*dy*sin(pi/2.d0-theta)/)
      end do

      if (layer.eq.2) then
        ver_hbox = ver_hbox2
      end if
      x_1 = ver_hbox(1,2)
      y_1 = ver_hbox(2,2)
      if (slope_bar .lt. 0.d0) then
        i_1 = i_0 - logical_to_int(xe(i_0).gt.x_1)
        j_1 = j_0 - logical_to_int(ye(j_0).gt.y_1)
      elseif (slope_bar .gt.0.d0) then
        i_1 = i_0 - logical_to_int(xe(i_0).gt.x_1) + logical_to_int(xe(i_0+1).lt.x_1)
        j_1 = j_0 - logical_to_int(ye(j_0).gt.y_1) + logical_to_int(ye(j_0+1).lt.y_1)
      end if


      print *, "Vers:",ver_hbox
      ! number of hboxes:
      m = size(intersect_top,2) - 1

      ! vertices of hboxes on the barrier
      allocate(we(2,size1))
      we = ver_hbox(:,odd)
      allocate(wc(2,m))
      do i =1,m
        wc(:,i) = (/0.5d0*(we(1,i+1)+we(1,i)),0.5d0*(we(2,i+1)+we(2,i))/)
      end do

      !bottom edges :
      allocate(we_1(2,size1))
      we_1 = ver_hbox(:,even)
      allocate(wc_1(2,m))
      do i=1,m
        wc_1(:,i)= (/0.5d0*(we_1(1,i+1)+we_1(1,i)),0.5d0*(we_1(2,i+1)+we_1(2,i))/)
      end do

      !hbox centers:
      allocate(hbox_cen(2,m))
      do i=1,m
        hbox_cen(:,i)=(/0.5d0*(wc(1,i)+wc_1(1,i)),0.5d0*(wc_1(2,i)+wc(2,i))/)
      end do

      ! allocating index array of small cells
      ! call cut_cells_inds(hbox_cen(1,:),hbox_cen_y(2,:),mx,my,mbc,xe,ye,ii,jj)

      ! for normal edges:
      if (slope_bar .lt. 0.d0 ) then
         allocate(i_chk(3))
         allocate(j_chk(3))
         i_chk = (/-1,0,1/)
         j_chk = (/-1,0,1/)
       elseif (slope_bar .gt. 0.d0) then
         allocate(i_chk(2))
         allocate(j_chk(2))
         i_chk = (/0,1/)
         j_chk = (/-1,0/)
       end if

       ! number of intersections in edges:
       allocate(no_int_norm(m+1))
       allocate(no_int_par(m))
       no_int_norm = 0
       no_int_par = 0
       if (slope_bar .lt. 0.d0) then
        in1 = 1
        in2 = 2
      else if ( slope_bar .gt. 0.d0) then
        in1 =  2
        in2 = 1
      end if
      allocate(i_range(size(i_chk,1)))
      allocate(j_range(size(j_chk,1)))
      i1 = 1 ! indexing counter for intersect_n
      do k=1,m+1
        xrange = (/ min(ver_hbox(1,2*k-1),ver_hbox(1,2*k)), max(ver_hbox(1,2*k-1),ver_hbox(1,2*k))/)
        yrange = (/min(ver_hbox(2,2*k-1),ver_hbox(2,2*k)),max(ver_hbox(2,2*k-1),ver_hbox(2,2*k))/)
        i_range = i_0_t + i_chk
        j_range = j_0_t + j_chk
        ci = 0
        do i=1,size(i_range,1)
          if (xe(i_range(i)).lt.xrange(2) .and. xe(i_range(i)).gt.xrange(1)) then
            dist_x = xe(i_range(i)) - xrange(1)
            dist_y = dist_x * tan(sign(1.d0,slope_bar)*(pi/2.d0-theta))
            dist_x2 = xe(i_range(i)) - xrange(2)
            if (abs(dist_x2).gt.tol .and. abs(dist_x).gt.tol) then
              call AddToList_verts(intersect_n,(/xe(i_range(i)),yrange(in1)-sign(1.d0,slope_bar)*dist_y/))
              no_int_norm(k) = no_int_norm(k) + 1
              ci = ci + 1
            end if
          end if
        end do
        do i = 1,size(j_range,1)
          if (ye(j_range(i)).lt.yrange(2) .and. ye(j_range(i)).gt.yrange(1)) then
            dist_y = abs(ye(j_range(i))-yrange(in1))
            dist_y2 = ye(j_range(i)) - yrange(in2)
            dist_x = dist_y / tan(sign(1.d0,slope_bar)*(pi/2.d0-theta))
            if (abs(dist_y2) .gt. tol .and. abs(dist_y) .gt. tol) then
              call AddToList_verts(intersect_n,(/xrange(1)+dist_x,ye(j_range(i))/))
              no_int_norm(k)= no_int_norm(k) + 1
              ci = ci+1
            end if
          end if
        end do

        if (ci .gt. 1) then
          do i = 1,ci
            call AddToList_verts(pair_inter_norm,intersect_n(:,size(intersect_n,2)-(i-1)))
          end do
        end if
        if (k.eq.m+1) THEN
          Exit
        endif
        if (slope_bar .lt. 0.d0 ) then
          i_0_t = i_0_t + logical_to_int(xe(i_0_t+1).lt.ver_hbox(1,2*k+1))
          j_0_t = j_0_t - logical_to_int(ye(j_0_t).gt.ver_hbox(2,2*k+1))
        else if (slope_bar .gt. 0.d0) then
          i_0_t = i_0_t + logical_to_int(xe(i_0_t+1).lt.ver_hbox(1,2*k+1))
          j_0_t = j_0_t + sign(1.d0,slope_bar)*logical_to_int(ye(j_0_t+1).lt.ver_hbox(2,2*k+1)) &
          -sign(1.d0,slope_bar) * logical_to_int(ye(j_0_t).gt.ver_hbox(2,2*k+1))
        end if
        if (k.eq.1) THEN
          bot_guess = (/i_0_t,j_0_t/)
        else if (k.eq.m) then
          bot_guess2= (/i_0_t,j_0_t/)
        end if
      end do

        ! for edges parallel to the barrier
        if (slope_bar.lt.0.d0) then
          allocate(j_chk2(3))
          j_chk2 = (/-1,0,1/)
          angle = pi - theta
        else if (slope_bar.gt.0.d0) then
          allocate(j_chk2(2))
          j_chk2 = (/0,1/)
          angle = theta
        end if
        allocate(i_range2(size(i_chk,1)))
        allocate(j_range2(size(j_chk2,1)))

        do k=1,m
          xrange = (/min(ver_hbox(1,2*k-1),ver_hbox(1,2*k+1)),max(ver_hbox(1,2*k-1),ver_hbox(1,2*k+1))/)
          xrange2 = (/min(ver_hbox(1,2*k),ver_hbox(1,2*k+2)),max(ver_hbox(1,2*k),ver_hbox(1,2*k+2))/)
          yrange =(/min(ver_hbox(2,2*k-1),ver_hbox(2,2*k+1)),max(ver_hbox(2,2*k-1),ver_hbox(2,2*k+1))/)
          yrange2 =(/min(ver_hbox(2,2*k),ver_hbox(2,2*k+2)),max(ver_hbox(2,2*k),ver_hbox(2,2*k+2))/)

          i_range = i_0_p + i_chk
          j_range = j_0_p + j_chk2

          ci=0
          do i=1,size(i_range,1)
            if (xe(i_range(i)).lt.xrange(2).and.xe(i_range(i)).gt.xrange(1)) then
              dist_x = abs(xrange(in2) - xe(i_range(i)))
              dist_x2 = xrange(in1) - xe(i_range(i))
              dist_y = dist_x * tan(angle)
              if (abs(dist_x2).gt.tol .and. abs(dist_x).gt.tol) then
                call AddToList_verts(intersect_p,(/xe(i_range(i)),yrange(1)+dist_y/))
                no_int_par(k) = no_int_par(k) + 1
                ci = ci + 1
              end if
            end if
          end do

          do i=1,size(j_range,1)
            if (ye(j_range(i)).lt.yrange(2) .and. ye(j_range(i)).gt.yrange(1)) then
              dist_y = abs(yrange(in2) - ye(j_range(i)))
              dist_y2 = yrange(in1) - ye(j_range(i))
              dist_x = dist_y/tan(angle)
              if (abs(dist_y2).gt.tol .and. abs(dist_y).gt.tol) then
                call AddToList_verts(intersect_p,(/xrange(1)+dist_x,ye(j_range(i))/))
                no_int_par(k) = no_int_par(k) + 1
                ci = ci + 1
              end if
            end if
          end do

          if (ci .gt. 1) then
            do i = 1,ci
              call AddToList_verts(pair_inter_par1,intersect_p(:,size(intersect_p,2)-(i-1)))
            end do
          end if
          if ( slope_bar.lt.0.d0) then
            i_0_p = i_0_p + logical_to_int(xe(i_0_p+1).lt.ver_hbox(1,2*k+1))
            j_0_p = j_0_p + sign(1.d0,slope_bar)*logical_to_int(ye(j_0_p+1).lt.ver_hbox(2,2*k+1)) &
            + sign(1.d0,slope_bar) * logical_to_int(ye(j_0_p).gt.ver_hbox(2,2*k+1))
          else if (slope_bar.gt.0.d0) then
            i_0_p = i_0_p + logical_to_int(xe(i_0_p+1).lt.ver_hbox(1,2*k+1))
            j_0_p = j_0_p + sign(1.d0,slope_bar)*logical_to_int(ye(j_0_p+1).lt.ver_hbox(2,2*k+1))&
            - sign(1.d0,slope_bar)*logical_to_int(ye(j_0_p).gt.ver_hbox(2,2*k+1))
          end if

          i_range2 = i_1 + i_chk
          j_range2 = j_1 + j_chk2

          ci=0
          do i = 1,size(i_range2,1)
            if (xe(i_range2(i)).lt.xrange2(2) .and. xe(i_range2(i)).gt.xrange2(1)) then
              dist_x2 = abs(xrange2(in2) - xe(i_range2(i)))
              dist_y2 = dist_x2 * tan(angle)
              dist_x3 = xrange2(in1) - xe(i_range2(i))
              if (abs(dist_x3).gt.tol .and. abs(dist_x2).gt.tol) then
                call AddToList_verts(intersect_p,(/xe(i_range2(i)),yrange2(1)+dist_y2/))
                no_int_par(k) = no_int_par(k) + 1
                ci = ci+1
              end if
            end if
          enddo

          do i =1,size(j_range2,1)
            if (ye(j_range2(i)).lt.yrange2(2) .and. ye(j_range2(i)).gt.yrange2(1)) then
              dist_y2 = abs(yrange2(in2) - ye(j_range2(i)))
              dist_x2 = dist_y2 / tan(angle)
              dist_y3 = yrange2(in1) - ye(j_range2(i))
              if (abs(dist_y3).gt.tol .and. abs(dist_y2).gt.tol) then
                call AddToList_verts(intersect_p,(/xrange2(1)+dist_x2,ye(j_range2(i))/))
                no_int_par(k) = no_int_par(k) + 1
                ci= ci+1
              end if
            end if
          end do

          if (ci .gt. 1) then
            do i = 1,ci
              call AddToList_verts(pair_inter_par2,intersect_p(:,size(intersect_p,2)-(i-1)))
            end do
          end if

          if (k.eq.m)then
            exit
          end if
          if (slope_bar.lt.0.d0) then
            i_1 = i_1 + logical_to_int(xe(i_1+1).lt.ver_hbox(1,2*k+2))
            j_1 = j_1 + sign(1.d0,slope_bar)*logical_to_int(ye(j_1+1).lt.ver_hbox(2,2*k+2)) &
            + sign(1.d0,slope_bar)*logical_to_int(ye(j_1).gt.ver_hbox(2,2*k+2))
          else if (slope_bar.gt.0.d0) then
            i_1 =  i_1 + logical_to_int(xe(i_1+1).lt.ver_hbox(1,2*k+2))
            j_1 = j_1 + sign(1.d0,slope_bar)*logical_to_int(ye(j_1+1).lt.ver_hbox(2,2*k+2))&
             - sign(1.d0,slope_bar)*logical_to_int(ye(j_1).gt.ver_hbox(2,2*k+2))
           end if
         end do

        ! physical vertices inside or not:
        i_h = i_0 - logical_to_int(xe(i_0).gt.hbox_cen(1,1)) &
         + logical_to_int(xe(i_0+1).gt.hbox_cen(1,1))
        j_h = j_0 - logical_to_int(ye(j_0).gt.hbox_cen(2,1)) + logical_to_int(ye(j_0+1).lt.hbox_cen(2,1))
        call IAddToList_verts(rough_ind,(/i_h,j_h/))

        allocate(phys_ver_in(m))
        ! print *, "m:",m
        ! print*, "pair inter p:", pair_inter_par1
        do k =1,m
          cands(1,:) = (/xe(i_h),xe(i_h),xe(i_h+1),xe(i_h+1)/)
          cands(2,:) = (/ye(j_h),ye(j_h+1),ye(j_h),ye(j_h+1)/)
          cen_x = hbox_cen(1,k)
          cen_y = hbox_cen(2,k)
          !hbox corners
          A = ver_hbox(:,2*k-1)
          B = ver_hbox(:,2*k)
          C = ver_hbox(:,2*k+1)
          D = ver_hbox(:,2*k+2)



          chk_in=0
          sides(1,:,:) = reshape((/A,B/),(/2,2/))
          sides(2,:,:) = reshape((/A,C/),(/2,2/))
          sides(3,:,:) = reshape((/C,D/),(/2,2/))
          sides(4,:,:) = reshape((/B,D/),(/2,2/))

          do i=1,4
            P = (/cands(1,i),cands(2,i)/)
            ! print*, "P:", P
            code = in_or_out(P,sides,4)
            ! print*," IN or not:", code
            if (code.eq.1 .and. dist(A,P).gt.tol .and. dist(B,P).gt.tol &
            .and. dist(C,P).gt.tol .and. dist(D,P).gt.tol) then
              chk_in = 1
              call AddToList_verts(phys_ver,P)
            end if
          end do
          if (chk_in.eq.1) then
            phys_ver_in(k) = 1
          else
            phys_ver_in(k) = 0
          end if
          if (k.eq.m) then
            exit
          end if
          i_h = i_h - logical_to_int(xe(i_h).gt.hbox_cen(1,k+1)) + &
            logical_to_int(xe(i_h+1).lt.hbox_cen(1,k+1))
          j_h = j_h + logical_to_int(ye(j_h+1).lt.hbox_cen(2,k+1)) -&
            logical_to_int(ye(j_h).gt.hbox_cen(2,k+1))
          call IAddToList_verts(rough_ind,(/i_h,j_h/))
        end do
        ! print*, "phys ver:", phys_ver

        ! Getting edges of hboxes
        if (allocated(pair_inter_norm)) then
          pair_inter_n = remove_dups(pair_inter_norm)
        endif
        if (allocated(pair_inter_par1)) then
          pair_inter_p1 = remove_dups(pair_inter_par1)
        end if
        if (allocated(pair_inter_par2)) then
          pair_inter_p2 = remove_dups(pair_inter_par2)
        end if

        allocate(no_int_norm2(m))
        no_int_norm2 = no_int_norm(2:) + no_int_norm(:m)
        ! index counter/pointers
        i = 0
        j = 0
        j2 = 0
        k1 = 0
        k2 = 0
        k3 = 0
        ! num of edges in hbox
        allocate(no_edges(m))
        ! addending edges
        do k = 1,m  !!!!!!!!!!!!!!!!!!!!!!!!!! ID 0
          ! print* , "K:", k
          num_edges = 0
          num_norm = no_int_norm2(k)
          num_par = no_int_par(k)
          A = ver_hbox(:,2*k-1)
          B = ver_hbox(:,2*k)
          C = ver_hbox(:,2*k+1)
          D = ver_hbox(:,2*k+2)
          x_vers = (/A(1),B(1),C(1),D(1)/)
          y_vers = (/A(2),B(2),C(2),D(2)/)
          call KB07AD(x_vers,4,dummy)
          call KB07AD(y_vers,4,dummy)
          xmin = x_vers(1)
          xmax = x_vers(4)
          ymin = y_vers(1)
          ymax = y_vers(4)

          normal_slope = ((C(2)-D(2))/(C(1)-D(1)))
          c_int = C(2) - normal_slope*C(1)
          do i2 = j,max(j+num_norm-1,0)
            if (i2.gt.size(intersect_n,2) .or. num_norm.eq.0) then
              exit
            end if
            if (allocated(intersect_n)) then
              call AddToList_verts(inter_chk,intersect_n(:,i2+1))
              call AddToList_verts(inter_chk_n,intersect_n(:,i2+1))
            end if
          end do

          do i2 = j2, max(j2+num_par-1,0)
            if (i2 .gt. size(intersect_p,2) .or. num_par.eq.0) then
              exit
            end if
            if (allocated(intersect_p)) then
              call AddToList_verts(inter_chk,intersect_p(:,i2+1))
            end if
          end do

          if (.not. allocated(inter_chk_n)) then
            j = j + num_norm
          else if (allocated(inter_chk_n)) then
            if (size(inter_chk_n,2).eq.1) then
              last_ver = inter_chk_n(:,1)
              if (abs(last_ver(2)-(normal_slope*last_ver(1) + c_int)).lt.tol) then
                j = j + num_norm - 1
              else
                j = j + 1
              end if
            end if
            if (size(inter_chk_n,2).gt.1) then
              last_ver = inter_chk_n(:,size(inter_chk_n,2))
              pen_ver = inter_chk_n(:,size(inter_chk_n,2)-1)
              if (abs(last_ver(2)-(normal_slope*last_ver(1)+c_int)).lt.tol &
              .and. abs(pen_ver(2)-(normal_slope*pen_ver(1) +c_int)).lt.tol) then
                j = j + num_norm - 2
              else if (abs(last_ver(2)-(normal_slope*last_ver(1)+c_int)).lt.tol) then
                j = j +num_norm -1
              else
                j = j +num_norm
              end if
            end if
          end if

          j2 = j2 + num_par

          ! remove duplicates:
          if (allocated(inter_chk)) then
            inter_chkd = remove_dups(inter_chk)
          end if
          print*,"intersections:", inter_chk

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EDGES SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! get all the right edges:
          if (.not. allocated(inter_chk)) then
            call AddToList_edges(edges,reshape((/A,B/),(/2,2/)))
            call AddToList_edges(edges,reshape((/A,C/),(/2,2/)))
            call AddToList_edges(edges,reshape((/C,D/),(/2,2/)))
            call AddToList_edges(edges,reshape((/B,D/),(/2,2/)))
            num_edges = num_edges + 4
          else
            allocate(wherea(size(inter_chkd,2),4))
            wherea = 0
            if (phys_ver_in(k).eq.1) then
              if (i.ge.size(phys_ver,2)) then
                i = i - 1
              end if
              x_phys_ver = phys_ver(1,i+1)
              y_phys_ver = phys_ver(2,i+1)
              ! print *, "x,y", x_phys_ver,y_phys_ver
              i = i + 1
              if (i.ge.size(phys_ver,2)) then
                i = i - 1
              end if
              call AddToList_verts(phizverts,(/x_phys_ver,y_phys_ver/))
              sides(1,:,:) = reshape((/A,B/),(/2,2/))
              sides(2,:,:) = reshape((/A,C/),(/2,2/))
              sides(3,:,:) = reshape((/C,D/),(/2,2/))
              sides(4,:,:) = reshape((/B,D/),(/2,2/))
              P = (/phys_ver(1,i+1),phys_ver(2,i+1)/)
              code = in_or_out(P,sides,4)
              if (code.eq.1) then
                call AddToList_verts(phizverts,(/phys_ver(1,i+1),phys_ver(2,i+1)/))
                i = i+1
              end if
                d_check = 1111111.d0
              if (allocated(phizverts)) then         !! need an end if statement somewhere down below
                ! if (size(phizverts,2).gt.1) then
                !   phiz2 = phizverts(:,2)
                !   phiz1 = phizverts(:,1)
                !   d_check = dist(phiz1,phiz2)
                !   if (d_check.lt.tol) then
                !     d_check = huge(0.d0)
                !     call DelFromList_verts(phizverts)
                !   end if
                !  end if
                 ! do p1=1,size(phizverts,2)  !! HERE 1
                    p1 = 1
                   x_phiz = phizverts(1,p1)
                   y_phiz = phizverts(2,p1)
                   do m1 = 1,size(inter_chkd,2) !! HERE 2
                      if (inter_chkd(1,m1) .gt. xmax .or. inter_chkd(1,m1).lt.xmin .or. &
                      inter_chkd(2,m1).gt. ymax .or. inter_chkd(2,m1).lt.ymin) then
                      continue
                elseif ((abs(x_phiz-inter_chkd(1,m1)).lt.tol.or.abs(y_phiz-inter_chkd(2,m1)).lt.tol)&
                      .and.(abs(x_phiz-inter_chkd(1,m1)).gt.tol.or.abs(y_phiz-inter_chkd(2,m1)).gt.tol)) then
                  if (size(phizverts,2).gt.1) then
                    if (dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).lt.d_check .and. &
                    dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).gt.tol) then
                      v1 = (/x_phiz,y_phiz/)
                      call AddToList_edges(edges,reshape((/v1,inter_chkd(:,m1)/),(/2,2/)))
                      num_edges = num_edges + 1
                    end if
                  else if (dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).gt.tol) then
                    v1 = (/x_phiz,y_phiz/)
                    call AddToList_edges(edges,reshape((/v1,inter_chkd(:,m1)/),(/2,2/)))
                    num_edges = num_edges + 1
                  endif
                end if
                normal_slope1 = (A(2)-B(2))/(A(1)-B(1))
                c1 = A(2) - normal_slope1*A(1)
                normal_slope2 = (A(2)-C(2))/(A(1)-C(1))
                c2 = A(2) - normal_slope2*A(1)
                normal_slope3 = (D(2)-B(2))/(D(1)-B(1))
                c3 = D(2) - normal_slope3*D(1)
                normal_slope4 = (D(2)-C(2))/(D(1)-C(1))
                c4 = D(2) - normal_slope4*D(1)
                if (abs(inter_chkd(2,m1)-(normal_slope1*inter_chkd(1,m1)+c1)).lt.tol) then
                  wherea(m1,1) = wherea(m1,1) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope2*inter_chkd(1,m1)+c2)).lt.tol) then
                  wherea(m1,2) = wherea(m1,2) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope3*inter_chkd(1,m1)+c3)).lt.tol) then
                  wherea(m1,3) = wherea(m1,3) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope4*inter_chkd(1,m1)+c4)).lt.tol) then
                  wherea(m1,4) = wherea(m1,4) + 1
                end if
              end do

              ! get which row has nonzero for each column, thats the index of the intersections.
              ! print *, "PHIZ:",x_phiz,y_phiz
              ! print *, " A:", A
              ! print *, " B: ",B
              ! print *, "C:", C
              ! print *, " D:",D
              if ((abs(x_phiz-A(1)).lt.tol .or. abs(y_phiz-A(2)).lt.tol) .and. &
                 dist(phizverts(:,p1),A).lt.d_check ) then
                 v1 = (/x_phiz,y_phiz/)
                 call AddToList_edges(edges,reshape((/v1,A/),(/2,2/)))
                 num_edges = num_edges + 1
               end if

               if ((abs(x_phiz-B(1)).lt.tol .or. abs(y_phiz-B(2)).lt.tol).and. &
                  dist(phizverts(:,p1),B).lt.d_check ) then
                  v1 = (/x_phiz,y_phiz/)
                  call AddToList_edges(edges,reshape((/v1,B/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
               if ((abs(x_phiz-C(1)).lt.tol .or. abs(y_phiz-C(2)).lt.tol) .and. &
                  dist(phizverts(:,p1),C).lt.d_check ) then
                  v1 = (/x_phiz,y_phiz/)
                  ! print *, " AM CALLED"
                  call AddToList_edges(edges,reshape((/v1,C/),(/2,2/)))
                  num_edges = num_edges + 1
                end if

               if ((abs(x_phiz-D(1)).lt.tol .or. abs(y_phiz-D(2)).lt.tol).and. &
                  dist(phizverts(:,p1),D).lt.d_check ) then
                  v1 = (/x_phiz,y_phiz/)
                  call AddToList_edges(edges,reshape((/v1,D/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
              ! end do  ! p1
            end if
               ! print*, "phiz vert", phizverts
               ! print *, "interchkd", inter_chkd
               ! print *, "size intchkd", size(inter_chkd,2)
              if (allocated(phizverts)) then  !#################   ID 1 need end if somewhere down
                if (size(phizverts,2).eq.1) then   !! ID 2
                  phiz = phizverts(:,1)
                  x_phiz = phiz(1)
                  y_phiz = phiz(2)
                  allocate(listt(size(inter_chkd,2)))
                  do m1 = 1,size(inter_chkd,2)

                    listt = (/(i2,i2=1,size(inter_chkd,2))/)
                    call IDelFromList(listt,m1) ! PROBELM IS IN THIS FUNCTION!!!
                    ! print *, "LIST", listt

                    do m2 = 1,size(listt,1)
                      m3 = listt(m2)
                      ! print *, "m3:",m3

                      if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                      abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                      (abs(inter_chkd(1,m1)-x_phiz).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz).gt.tol)&
                      .and.(abs(inter_chkd(1,m1)-inter_chkd(1,m3)).gt.tol.or.abs(inter_chkd(2,m1)-&
                      inter_chkd(2,m3)).gt.tol)) then
                        call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                        num_edges = num_edges + 1
                      end if
                    end do

                  end do

                ! else if (size(phizverts,2).gt.1) then   !! ID 2 IF CLAUSE
                !   phiz1 = phizverts(:,1)
                !   x_phiz1 = phiz1(1)
                !   y_phiz1 = phiz1(2)
                !   phiz2 = phizverts(:,2)
                !   x_phiz2 = phiz2(1)
                !   y_phiz2 = phiz2(2)
                !   do m1 = 1, size(inter_chkd,2)
                !     listt = (/(i2,i2=1,size(inter_chkd,2))/)
                !     call IDelFromList(listt,m1)
                !     do m2 = 1,size(listt,1)
                !       m3 = listt(m2)
                !       if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                !       abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                !       (abs(inter_chkd(1,m1)-x_phiz1).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz1).gt.tol)&
                !  .and.(abs(inter_chkd(1,m1)-x_phiz2).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz2).gt.tol)&
                !       .and.(abs(inter_chkd(1,m1)-inter_chkd(1,m3)).gt.tol.or.abs(inter_chkd(2,m1)-&
                !       inter_chkd(2,m3)).gt.tol)) then
                !       call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                !       num_edges = num_edges + 1
                !       end if
                !     end do
                !   end do
                !   if ((abs(x_phiz1-x_phiz2).lt.tol.or.abs(y_phiz1-y_phiz2).lt.tol).and. &
                !      dist(phiz1,phiz2).gt.tol) then
                !      m3 = size(phizverts,2)
                !      call AddToList_edges(edges,reshape((/phizverts(:,1),phizverts(:,m3)/),(/2,2/)))
                !      num_edges = num_edges + 1
                !   end if
                end if
              end if


            else  !!!!  ID 3
              ! print *, "here! no phiz vert"
              do m1 = 1,size(inter_chkd,2)
                listt = (/(i2,i2=1,size(inter_chkd,2))/)
                call IDelFromList(listt,m1)
                do m2 = 1,size(listt,1)
                  m3 = listt(m2)
                  if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                  abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                     dist(inter_chkd(:,m1),inter_chkd(:,m3)).gt.tol) then
                    call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                end do
                normal_slope1 = (A(2)-B(2))/(A(1)-B(1))
                c1 = A(2) - normal_slope1*A(1)
                normal_slope2 = (A(2)-C(2))/(A(1)-C(1))
                c2 = A(2) - normal_slope2*A(1)
                normal_slope3 = (D(2)-B(2))/(D(1)-B(1))
                c3 = D(2) - normal_slope3*D(1)
                normal_slope4 = (D(2)-C(2))/(D(1)-C(1))
                c4 = D(2) - normal_slope4*D(1)
                if (abs(inter_chkd(2,m1)-(normal_slope1*inter_chkd(1,m1)+c1)).lt.tol) then
                  wherea(m1,1) = wherea(m1,1) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope2*inter_chkd(1,m1)+c2)).lt.tol) then
                  wherea(m1,2) = wherea(m1,2) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope3*inter_chkd(1,m1)+c3)).lt.tol) then
                  wherea(m1,3) = wherea(m1,3) + 1
                else if (abs(inter_chkd(2,m1)-(normal_slope4*inter_chkd(1,m1)+c4)).lt.tol) then
                  wherea(m1,4) = wherea(m1,4) + 1
                end if
              end do
            end if

            x_cat(1,:) = (/A(1),A(1),B(1),C(1)/)
            x_cat(2,:) = (/B(1),C(1),D(1),D(1)/)
            y_cat(1,:) = (/A(2),A(2),B(2),C(2)/)
            y_cat(2,:) = (/B(2),C(2),D(2),D(2)/)
            where_sum = sum(wherea,dim=1)
            ! print *, "corners",A,B,C,D
            ! print *, "inter chdk",inter_chkd
            ! print*, "where", where_sum
            do m3 = 1,4
              v1 = (/x_cat(1,m3),y_cat(1,m3)/)
              v2 = (/x_cat(2,m3),y_cat(2,m3)/)
              if ( where_sum(m3).eq.0 .and.dist(v1,v2).gt.tol) then
               call AddToList_edges(edges,reshape((/v1,v2/),(/2,2/)))
               num_edges = num_edges + 1
              end if
              if (where_sum(m3) .eq. 1) then
                o = nonzero(wherea(:,m3))
                if (dist(inter_chkd(:,o),v1).gt.tol)THEN
                  call AddToList_edges(edges,reshape((/inter_chkd(:,o),v1/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
                if (dist(inter_chkd(:,o),v2).gt.tol)then
                  call AddToList_edges(edges,reshape((/inter_chkd(:,o),v2/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
              end if
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1333
              if (where_sum(m3) .gt. 1) then
                find = nonzero(wherea(:,m3))
                ! print *, "fnd!",find
                allocate(x_list(size(find)))
                x_list = inter_chkd(1,find)
                min_loc = minloc(x_list)
                ! print *, "xlist",x_list
                ml = min_loc(1)
                ! print*,"Min_loc:",  ml
                max_loc = maxloc(x_list)
                nl = max_loc(1)
                ver_x = x_cat(:,m3)
                ver_min = minloc(ver_x)
                vl = ver_min(1)
                ver_max = maxloc(ver_x)
                wl = ver_max(1)
                v3 = inter_chkd(:,find(ml))
                v4 = (/minval(x_cat(:,m3)),y_cat(vl,m3)/)
                ! print*, "dbl eddg1", v3,v4
                if (dist(v3,v4).gt.tol) then
                  call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
                v3 = inter_chkd(:,find(nl))
                v4 = (/maxval(x_cat(:,m3)),y_cat(wl,m3)/)
                ! print*, "dbl eddg2", v3,v4
                if (dist(v3,v4).gt.tol) then
                  call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
                v4 = inter_chkd(:,find(ml))
                ! print*, "dbl eddg3", v3,v4
                if (dist(v3,v4).gt.tol) then
                  call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                  num_edges = num_edges + 1
                end if
                deallocate(x_list)
              end if
            end do
          end if

          no_edges(k) =  num_edges
          if (allocated(wherea)) then
            deallocate(wherea)
          end if
          if (allocated(phizverts)) then
            deallocate(phizverts)
          end if
          if (allocated(listt)) then
            deallocate(listt)
          end if
          if (allocated(inter_chkd)) then
            deallocate(inter_chkd)
            deallocate(inter_chk)
          end if
        end do

        !!!!!!!!!!!1! getting the fragments of each hbox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        counter = 0
        counter2 = 0
        ! allocate(q_hbox(3,m))
        ! allocate(aux_hbox(1,m))
        ! q_hbox = 0.d0
        ! aux_hbox = 0.d0
        ! print*, "edgesL",edges
        do k = 1, m
          print*, " K: ",k
          num_edge = no_edges(k)

          !!!!!!!!!!!*************PROBLEM IS DOWN HERE WITH THE RIGHT INDEX FINDING FUNCTIONSSS**************!!!!!!!!!!!!!!!!!

          do j=counter,max(counter+num_edge-1,0)
            e1 = edges(:,:,j+1)
            ord = vec2ord(e1(1,:))
            call AddToList_edges(edges_chk,e1(:,ord))
          end do
          ! deallocate(indexes)
          ! allocate(indexes(size(edges_chk,3)))
          ! call KB07AD(edges_chk(1,1,:),size(edges_chk,3),indexes)
          ! edges_chk(:,:,1:size(edges_chk,3)) = edges_chk(:,:,indexes)
          counter = counter + num_edge
          print*, "Edges: ordered ", edges_chk
          ! hbox area info and list to collect areas and indices of frags
          ! hbox vertices
          A = ver_hbox(:,2*k-1)
          B = ver_hbox(:,2*k)
          C = ver_hbox(:,2*k+1)
          D = ver_hbox(:,2*k+2)
          hbox_area = area_polygon((/A(1),B(1),D(1),C(1),A(1)/),(/A(2),B(2),D(2),C(2),A(2)/) )!dist(A,B)*dist(B,D)
          print *, "HBOX AREA: ",hbox_area

          init_area = 0.d0

          if (num_edge .eq. 4) then
            call AddToList_verts(cyclea,A)
            call AddToList_verts(cyclea,B)
            call AddToList_verts(cyclea,D)
            call AddToList_verts(cyclea,C)
            call AddToList_verts(cyclea,A)
            call AddToList(areas,hbox_area)
            call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
            ,xe,ye,mx,my,mbc,right_ind,centroid)
            call IAddToList_verts(inds,right_ind)
            init_area = init_area + hbox_area
            num_frags(k) = 1
            hbox_areas(k) = hbox_area
            area_frags(1,k) = 1.d0
            index_frags(:,1,k) = right_ind
            goto 99
          end if
          if (allocated(cyclea)) then
            deallocate(cyclea)
          end if
          ! get the lateral two edges, "start edges" to do the rotation
          do j = 1, size(edges_chk,3)
            e1 = edges_chk(:,:,j)
            v4 = e1(:,2)
            v5 = e1(:,1)
            v6 = v4 - v5
            if (abs(v6(2)) .lt. tol) then
              if (.not. allocated(horz_edges)) then
                call AddToList_edges(horz_edges,e1)
              else ! check if edge e1 is not in horz edges and e1 permuted not in horz edges
                if (.not. in(e1,horz_edges)) then
                  call AddToList_edges(horz_edges,e1)
                end if
              end if
            end if
            if (abs(v6(1)).lt.tol) then
              if (.not. allocated(vert_edges)) then
                call AddToList_edges(vert_edges,e1)
              else
                if (.not. in(e1,vert_edges)) then
                  call AddToList_edges(vert_edges,e1)
                end if
              end if
            end if
          end do
          ! if both horz and vert edges exist
          if (allocated(horz_edges)) then
            print*, "HORZ EDGES!",horz_edges
            ! horz edge counter/clkwise cycle search
            isign = -1
            do j=1,size(horz_edges,3)
              h_edge = horz_edges(:,:,j)
              do i = 1,2
                call AddToList_edges(region_d,h_edge)
                call AddToList_verts(cyclea,h_edge(:,1))
                call AddToList_verts(cyclea,h_edge(:,2))
                !!!!!!!!!!! here 9.18.20 5:38 PM
                call is_cycle(reshape(cyclea,(/2,2,1/)),yes_or_no)
                i2 = 1
                do while (.not.yes_or_no)
                  do i3=1,size(edges_chk,3)
                    e1 = edges_chk(:,:,i3)
                    v1 = e1(:,1)
                    v2 = e1(:,2)
                    x1 = v1(1)
                    x2 = v2(1)
                    y1 = v1(2)
                    y2 = v2(2)
                    last_elt = region_d(:,:,size(region_d,3))
                    v_last1 = last_elt(:,1)
                    v_last2 = last_elt(:,2)

                    if ((dist(v1,v_last2).lt.tol .or. dist(v2,v_last2).lt.tol).or.&
                    (dist(v1,v_last1).lt.tol .or. dist(v2,v_last1).lt.tol)) then
                       if (dist(v1,v_last2).lt.tol) then
                          cand = v2
                        else if (dist(v2,v_last2).lt.tol) then
                          cand = v1
                        else if (dist(v1,v_last1).lt.tol) then
                          cand = v2
                        else if (dist(v2,v_last1).lt.tol) then
                          cand = v1
                        end if

                        ver3 = cand
                        ver2 = cyclea(:,i2+1)
                        ver1 = cyclea(:,i2)
                        vec1 = ver3 - ver2
                        vec2 = ver2 - ver1

                        z = vec2(1) * vec1(2) - vec1(1)*vec2(2)
                        if (sign(1.d0,z*isign).gt.0.d0 .and. abs(z) .gt.tol**2 ) then
                          if (.not. in(e1, region_d).and..not.vert_in_edge(cyclea(:,i2),e1))&
                          then
                            call AddToList_verts(cyclea,cand)
                            call AddToList_edges(region_d,e1)
                            call is_cycle(region_d,yes_or_no)
                            i2 = i2 + 1
                            exit
                          end if
                        endif
                        continue
                      end if
                    end do
                    call is_cycle(region_d,yes_or_no)
                  end do ! while
                  area = area_polygon(cyclea(1,:),cyclea(2,:))
                  print *, " AREA" , area
                  print* , cyclea
                  call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
                   ,xe,ye,mx,my,mbc,right_ind,centroid)
                   print *, "right ind:",right_ind
                  if (.not. allocated(areas)) then
                    init_area = init_area + abs(area)
                    call AddToList(areas,area)
                    call IAddToList_verts(inds,right_ind)
                  else if (.not. ind_in_inds(right_ind,inds)) THEN
                    init_area= init_area + area
                    call AddToList(areas,area)
                    call IAddToList_verts(inds,right_ind)
                    if (abs(init_area-hbox_area).lt.tol*100.d0) then
                      goto 99
                    end if
                  end if
                  isign = isign * -1
                  deallocate(region_d)
                  deallocate(cyclea)
                end do
              end do
            end if
            print *, "Inds:", inds
            ! VERT EDGES :
            if (allocated(vert_edges)) then
              print*, "VERT EDGES!2", vert_edges
              ! horz edge counter/clkwise cycle search
              isign = -1
              do j=1,size(vert_edges,3)
                v_edge = vert_edges(:,:,j)
                do i = 1,2
                  call AddToList_edges(region_d,v_edge)
                  call AddToList_verts(cyclea,v_edge(:,1))
                  call AddToList_verts(cyclea,v_edge(:,2))
                  !!!!!!!!!!! here 9.18.20 5:38 PM
                  call is_cycle(reshape(cyclea,(/2,2,1/)),yes_or_no)
                  i2 = 1
                  do while (.not.yes_or_no)
                    do i3=1,size(edges_chk,3)
                      e1 = edges_chk(:,:,i3)
                      v1 = e1(:,1)
                      v2 = e1(:,2)
                      x1 = v1(1)
                      x2 = v2(1)
                      y1 = v1(2)
                      y2 = v2(2)
                      last_elt = region_d(:,:,size(region_d,3))
                      v_last1 = last_elt(:,1)
                      v_last2 = last_elt(:,2)

                      if ((dist(v1,v_last2).lt.tol .or. dist(v2,v_last2).lt.tol).or.&
                      (dist(v1,v_last1).lt.tol .or. dist(v2,v_last1).lt.tol)) then
                         if (dist(v1,v_last2).lt.tol) then
                            cand = v2
                          else if (dist(v2,v_last2).lt.tol) then
                            cand = v1
                          else if (dist(v1,v_last1).lt.tol) then
                            cand = v2
                          else if (dist(v2,v_last1).lt.tol) then
                            cand = v1
                          end if

                          ver3 = cand
                          ver2 = cyclea(:,i2+1)
                          ver1 = cyclea(:,i2)
                          vec1 = ver3 - ver2
                          vec2 = ver2 - ver1

                          z = vec2(1) * vec1(2) - vec1(1)*vec2(2)
                          if (sign(1.d0,z*isign).gt.0.d0 .and. abs(z) .gt.tol**2) then
                            if (.not. in(e1, region_d).and..not.vert_in_edge(cyclea(:,i2),e1))&
                            then
                              call AddToList_verts(cyclea,cand)
                              call AddToList_edges(region_d,e1)
                              call is_cycle(region_d,yes_or_no)
                              i2 = i2 + 1
                              exit
                            end if
                          endif
                          continue
                        end if
                      end do
                    end do ! while
                    area = area_polygon(cyclea(1,:),cyclea(2,:))
                    print *, " AREA" , area
                    print*, cyclea
                    call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
                     ,xe,ye,mx,my,mbc,right_ind,centroid)
                     print *, "RI",right_ind
                    if (.not. allocated(areas)) then
                      init_area = init_area + abs(area)
                      call AddToList(areas,area)
                      call IAddToList_verts(inds,right_ind)
                    else if (.not. ind_in_inds(right_ind,inds)) THEN
                      init_area= init_area + area
                      call AddToList(areas,area)
                      call IAddToList_verts(inds,right_ind)
                      if (abs(init_area-hbox_area).lt.tol*100.d0) then
                        goto 99
                      end if
                    end if
                    isign = isign * -1
                    deallocate(region_d)
                    deallocate(cyclea)
                  end do
                end do
              end if
              ! print*, "ri:",right_ind
     99 continue
              print *, "Areas: ", areas, hbox_area
              num_frags(k) = size(areas,1)
              hbox_areas(k) = hbox_area
              if (abs(init_area-hbox_area).gt.tol) THEN
                print *, "your code didn't work!!"
              end if
              print*," INDS", inds
              do i=1,size(areas,1)
                 look = inds(:,i)
                 area_frags(i,k) = areas(i)/hbox_area
                 index_frags(:,i,k) = look
                 ! q_hbox(:,k) = q_hbox(:,k) + q_old(:,look(1),look(2))*areas(i)/hbox_area
                 ! aux_hbox(k) = aux_hbox(k) + aux(1,look(1),look(2))*areas(i)/hbox_area
               enddo
               deallocate(edges_chk)
               if (allocated(horz_edges)) then
                 deallocate(horz_edges)
               end if
               if (allocated(vert_edges)) then
                 deallocate(vert_edges)
               end if
               if (allocated(areas)) then
                 deallocate(areas)
               end if
               if (allocated(inds)) then
                 deallocate(inds)
               endif
               if (allocated(region_d)) then
                 deallocate(region_d)
               end if
               if (allocated(cyclea)) then
                 deallocate(cyclea)
               end if
             end do

  end subroutine
    function logical_to_int(L)
      implicit none
      integer:: logical_to_int
      logical :: L
      if (L) then
        logical_to_int = 1
      else
        logical_to_int = 0
      end if
    end function

    subroutine up_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,i_0,j_0,x_0,y_0,i_e,j_e,x_e,y_e,&
      m,hbox_areas,num_frags,area_frags,index_frags,layer)!,ii,jj)

      implicit none
      ! input:
      real(8):: xlower,ylower,xupper,yupper,x_0,y_0,x_e,y_e ! xlower,xupper,ylower,yupper is the grid boundary points; x_0,y_0,x_e,y_e endpoints of barrier segment
      integer :: mbc,mx,my,i_0,j_0,i_e,j_e,layer !mbc,mx,my are num ghost, grid cells; i_0,j_0,i_e,j_e indices of endpoints of barrier
                                                 !layer is integer indicating first or second layer of hboxes
      real(8) :: q_old(3,1-mbc:mx+mbc,1-mbc:my+mbc) ! q_old is q values on the grid at time t^n, but this should be q_old_upper as opposed to q_old_lower
                                                    ! q_old_lower is the grid q values where the cut cells are replaced with the bottom small cell values
                                                    ! q_old_upper is the grid q values where the cut cells are replaced with the top small cell values
                                                    ! this is to simplify getting the hbox averages by indexing in the context of grid
      real(8) :: aux(1,1-mbc:mx+mbc,1-mbc:my+mbc)

      ! output:
      integer :: m  ! m is number of hboxes
      integer :: num_frags(mx*4), index_frags(2,4,mx*4) ! number of fragments in each hbox ; the index location of each fragment (<=4)
      real(8):: hbox_areas(mx*4), area_frags(4,mx*4)  ! area of each hbox ; areas (<=4) of ea frag in ea hbox

      ! local:
      real(8) :: dx, dy, xe(-mbc:mx+mbc), ye(-mbc:my+mbc),x1,x2,y1,y2
      real(8) :: slope_bar,theta,dist_x,dist_y,dist_x2,dist_y2,dist_x3
      real(8) :: x_1,y_1, xrange(2),yrange(2),angle,xrange2(2),yrange2(2)
      real(8) :: cands(2,4),cen_x,cen_y,A(2),B(2),C(2),D(2),sides(4,2,2),P(2)
      real(8) :: x_vers(4), y_vers(4),xmin,xmax,ymin,ymax,normal_slope, c_int
      real(8) :: last_ver(2), pen_ver(2),x_phys_ver,y_phys_ver,phiz1(2),phiz(2)
      real(8) :: d_check,x_phiz,y_phiz,v1(2),v2(2),normal_slope1,normal_slope2
      real(8) :: normal_slope3,normal_slope4,c1,c2,c3,c4,phiz2(2),x_phiz1,x_phiz2
      real(8) :: y_phiz1,y_phiz2,x_cat(2,4),y_cat(2,4),ver_x(2),v3(2),v4(2),v5(2),v6(2)
      real(8) :: e1(2,2),hbox_area,init_area,centroid(2),h_edge(2,2),last_elt(2,2),z
      real(8) :: v_last1(2),v_last2(2),cand(2),ver3(2),ver2(2),ver1(2),vec1(2),vec2(2)
      real(8) :: area,v_edge(2,2),dist_y3
      real(8),allocatable ::intersect_0(:,:), pair_inter_norm(:,:),pair_inter_par1(:,:)
      real(8),allocatable :: intersect_n(:,:),intersect_p(:,:),pair_inter_par2(:,:)
      real(8),allocatable :: intersect_top(:,:), ver_hbox(:,:),ver_hbox2(:,:),phizverts(:,:)
      real(8),allocatable :: we(:,:), wc(:,:), we_1(:,:), wc_1(:,:),hbox_cen(:,:),x_list(:)
      real(8),allocatable :: pair_inter_n(:,:), pair_inter_p1(:,:),pair_inter_p2(:,:)
      real(8),allocatable :: edges(:,:,:),inter_chk_n(:,:),inter_chk(:,:),inter_chkd(:,:)
      real(8),allocatable :: edges_chk(:,:,:),areas(:), cyclea(:,:),horz_edges(:,:,:)
      real(8),allocatable :: vert_edges(:,:,:),region_d(:,:,:),phys_ver(:,:)

      integer :: i,i2,j,j2,k1,k2,k3,i_0_p,i_0_t,j_0_p,j_0_t,k,ci,i0,r0,size1,num_edges
      integer :: num_norm, num_par,p1,m1,m2,m3,ml,nl,vl,wl
      integer:: i_1,j_1,in1,in2,i1,bot_guess(2),bot_guess2(2),i_h,j_h,chk_in,code
      integer :: counter,counter2,num_edge,ord(2), right_ind(2),isign,i3,look(2)
      integer :: slope_sign,i_chk_2(4),where_sum(4),dummy(4)
      integer,allocatable :: rough_ind(:,:),wherea(:,:),listt(:),find(:),o(:)
      integer , allocatable :: indexes(:), odd(:),even(:),i_chk(:),j_chk(:),i_range(:),j_range(:)
      integer, allocatable :: no_int_norm(:), no_int_par(:),i_chk2(:),j_chk2(:),i_range2(:),j_range2(:)
      integer , allocatable :: no_int_norm2(:),no_edges(:),inds(:,:)!, ii(:),jj(:)
      integer , allocatable :: phys_ver_in(:),min_loc(:),max_loc(:),ver_min(:),ver_max(:)

      logical :: vertical,yes_or_no

      area_frags=0.d0
      hbox_areas = 0.d0

      dx = (xupper-xlower)/mx
      dy = (yupper-ylower)/my
      xe = (/(xlower+i*dx, i=-2,mx+2)/)
      ye = (/(ylower+i*dy,i=-2,my+2)/)

      i_0_p = i_0
      i_0_t = i_0
      j_0_p = j_0
      j_0_t = j_0

      ! wall parameters:
      vertical = .false.
      if (x_e-x_0 .ne. 0.d0) then
        slope_bar = (y_e-y_0)/(x_e-x_0)
        theta = atan(slope_bar)
      else
        vertical = .true.
      end if
      slope_sign = sign(1.d0,slope_bar)

      ! get intersections between barrier and grid lines :
      do k=-2,mx+2
        if (xe(k).le. x_e+tol .and. xe(k).ge. x_0-tol) then
          dist_x = x_e - xe(k)
          dist_y = dist_x * tan(pi-theta)
          call AddToList_verts(intersect_0,(/xe(k),y_e+dist_y/))
        end if
      end do

      do k=-2,my+2
        if (ye(k).le.max(y_e,y_0)+tol .and. ye(k).ge.min(y_0,y_e)-tol) then
          dist_y = y_0 - ye(k)
          dist_x = dist_y/tan(pi-theta)
          call AddToList_verts(intersect_0,(/x_0+dist_x,ye(k)/))
        end if
      end do

      ! treat array (get rd of duplicate and sort)
         intersect_top = remove_dups(intersect_0)
         ! print*, "Intersections: x", intersect_top(1,:)
         ! print*, "intersections: y", intersect_top(2,:)
         size1 = size(intersect_top,2)
        allocate(indexes(size1))
        call KB07AD(intersect_top(1,:),size1,indexes)
        intersect_top(2,:) = intersect_top(2,indexes)
        ! print*, "Intersections: x", intersect_top(1,:)
        ! print*, "intersections: y", intersect_top(2,:)

        ! first layer of hbox vertices
        allocate(ver_hbox(2,2*size1))
        allocate(ver_hbox2(2,2*size1))

        allocate(odd(size1))
        allocate(even(size1))
        odd = (/ (i,i=1,2*size1,2)/)
        even = (/(i,i=2,2*size1,2)/)
        ver_hbox(:,odd) = intersect_top
        do i = 1,size1
          ver_hbox(:,even(i)) = intersect_top(:,i) + &
          (/-dx*cos(slope_sign*(pi/2.d0-theta)),slope_sign*dy*sin(slope_sign*(pi/2.d0-theta))/)
        end do

        ver_hbox2(:,odd) = ver_hbox(:,even)
        do i= 1,size1
          ver_hbox2(:,even(i)) = intersect_top(:,i) + &
          (/-2*dx*cos(slope_sign*(pi/2.d0-theta)),2*slope_sign*dy*sin(slope_sign*(pi/2.d0-theta))/)
        end do

        if (layer.eq.2) then
          ver_hbox = ver_hbox2
        end if

        x_1 = ver_hbox(1,2)
        y_1 = ver_hbox(2,2)
        if (slope_bar .lt. 0.d0) then
          i_1 = i_0 + logical_to_int(xe(i_0+1).lt.x_1)+1
          j_1 = j_0 + logical_to_int(ye(j_0+1).lt.y_1)+1
        elseif (slope_bar .gt.0.d0) then
          i_1 = i_0 - logical_to_int(xe(i_0).gt.x_1)+1
          j_1 = j_0 + logical_to_int(ye(j_0+1).lt.y_1)+1
        end if

        print*,"vers:", ver_hbox
        ! number of hboxes:
        m = size(intersect_top,2) - 1

        ! vertices of hboxes on the barrier
        allocate(we(2,size1))
        we = ver_hbox(:,odd)
        allocate(wc(2,m))
        do i =1,m
          wc(:,i) = (/0.5d0*(we(1,i+1)+we(1,i)),0.5d0*(we(2,i+1)+we(2,i))/)
        end do

        !bottom edges :
        allocate(we_1(2,size1))
        we_1 = ver_hbox(:,even)
        allocate(wc_1(2,m))
        do i=1,m
          wc_1(:,i)= (/0.5d0*(we_1(1,i+1)+we_1(1,i)),0.5d0*(we_1(2,i+1)+we_1(2,i))/)
        end do

        !hbox centers:
        allocate(hbox_cen(2,m))
        do i=1,m
          hbox_cen(:,i)=(/0.5d0*(wc(1,i)+wc_1(1,i)),0.5d0*(wc_1(2,i)+wc(2,i))/)
        end do

        ! allocating index array of small cells
        ! call cut_cells_inds(hbox_cen(1,:),hbox_cen_y(2,:),mx,my,mbc,xe,ye,ii,jj)

        ! for normal edges:
        if (slope_bar .lt. 0.d0 ) then
           allocate(i_chk(3))
           allocate(j_chk(3))
           i_chk = (/0,1,2/)
           j_chk = (/0,1,2/)
         elseif (slope_bar .gt. 0.d0) then
           allocate(i_chk(2))
           allocate(j_chk(2))
           i_chk = (/-1,0/)
           j_chk = (/0,1/)
         end if

         ! number of intersections in edges:
         allocate(no_int_norm(m+1))
         allocate(no_int_par(m))
         no_int_norm = 0
         no_int_par = 0
         if (slope_bar .lt. 0.d0) then
          in1 = 1
          in2 = 2
        else if ( slope_bar .gt. 0.d0) then
          in1 =  2
          in2 = 1
        end if
        allocate(i_range(size(i_chk,1)))
        allocate(j_range(size(j_chk,1)))
        i1 = 1 ! indexing counter for intersect_n
        do k=1,m+1
          xrange = (/ min(ver_hbox(1,2*k-1),ver_hbox(1,2*k)), max(ver_hbox(1,2*k-1),ver_hbox(1,2*k))/)
          yrange = (/min(ver_hbox(2,2*k-1),ver_hbox(2,2*k)),max(ver_hbox(2,2*k-1),ver_hbox(2,2*k))/)
          i_range = i_0_t + i_chk
          j_range = j_0_t + j_chk
          ci = 0
          do i=1,size(i_range,1)
            if (xe(i_range(i)).lt.xrange(2) .and. xe(i_range(i)).gt.xrange(1)) then
              dist_x = xe(i_range(i)) - xrange(1)
              dist_y = dist_x * tan(sign(1.d0,slope_bar)*(pi/2.d0-theta))
              dist_x2 = xe(i_range(i)) - xrange(2)
              if (abs(dist_x2).gt.tol .and. abs(dist_x).gt.tol) then
                call AddToList_verts(intersect_n,(/xe(i_range(i)),yrange(in1)-sign(1.d0,slope_bar)*dist_y/))
                no_int_norm(k) = no_int_norm(k) + 1
                ci = ci + 1
              end if
            end if
          end do
          do i = 1,size(j_range,1)
            if (ye(j_range(i)).lt.yrange(2) .and. ye(j_range(i)).gt.yrange(1)) then
              dist_y = abs(ye(j_range(i))-yrange(in1))
              dist_y2 = ye(j_range(i)) - yrange(in2)
              dist_x = dist_y / tan(sign(1.d0,slope_bar)*(pi/2.d0-theta))
              if (abs(dist_y2) .gt. tol .and. abs(dist_y) .gt. tol) then
                call AddToList_verts(intersect_n,(/xrange(1)+dist_x,ye(j_range(i))/))
                no_int_norm(k)= no_int_norm(k) + 1
                ci = ci+1
              end if
            end if
          end do

          if (ci .gt. 1) then
            do i = 1,ci
              call AddToList_verts(pair_inter_norm,intersect_n(:,size(intersect_n,2)-(i-1)))
            end do
          end if
          if (k.eq.m+1) THEN
            Exit
          endif
          if (slope_bar .lt. 0.d0 ) then
            i_0_t = i_0_t + logical_to_int(xe(i_0_t+1).lt.ver_hbox(1,2*k+1))
            j_0_t = j_0_t - logical_to_int(ye(j_0_t).gt.ver_hbox(2,2*k+1))
          else if (slope_bar .gt. 0.d0) then
            i_0_t = i_0_t + logical_to_int(xe(i_0_t+1).lt.ver_hbox(1,2*k+1)) &
                - logical_to_int(xe(i_0_t).gt.ver_hbox(1,2*k+1))
            j_0_t = j_0_t + sign(1.d0,slope_bar)*logical_to_int(ye(j_0_t+1).lt.ver_hbox(2,2*k+1)) &
            -sign(1.d0,slope_bar) * logical_to_int(ye(j_0_t).gt.ver_hbox(2,2*k+1))
          end if
          if (k.eq.1) THEN
            bot_guess = (/i_0_t,j_0_t/)
          else if (k.eq.m) then
            bot_guess2= (/i_0_t,j_0_t/)
          end if
        end do

          ! for edges parallel to the barrier
          i_chk_2 = (/-1,0,1,2/)
          if (slope_bar.lt.0.d0) then
            allocate(j_chk2(2))
            j_chk2 = (/-1,0/)
            angle = pi - theta
          else if (slope_bar.gt.0.d0) then
            allocate(j_chk2(3))
            j_chk2 = (/0,1,2/)
            angle = theta
          end if
          deallocate(i_range)
          deallocate(j_range)
          allocate(i_range(size(i_chk_2)))
          allocate(j_range(size(j_chk2)))
          allocate(i_range2(4))
          allocate(j_range2(size(j_chk2,1)))

          do k=1,m
            xrange = (/min(ver_hbox(1,2*k-1),ver_hbox(1,2*k+1)),max(ver_hbox(1,2*k-1),ver_hbox(1,2*k+1))/)
            xrange2 = (/min(ver_hbox(1,2*k),ver_hbox(1,2*k+2)),max(ver_hbox(1,2*k),ver_hbox(1,2*k+2))/)
            yrange =(/min(ver_hbox(2,2*k-1),ver_hbox(2,2*k+1)),max(ver_hbox(2,2*k-1),ver_hbox(2,2*k+1))/)
            yrange2 =(/min(ver_hbox(2,2*k),ver_hbox(2,2*k+2)),max(ver_hbox(2,2*k),ver_hbox(2,2*k+2))/)

            i_range = i_0_p + i_chk_2
            j_range = j_0_p + j_chk2

            ci=0
            print *, "k: ",k
            print*,"range: x",xrange
            print*,"range : y",yrange

            print*,"range: x2",xrange2
            print*,"range : y2",yrange2
            do i=1,size(i_range,1)

              if (xe(i_range(i)).lt.xrange(2).and.xe(i_range(i)).gt.xrange(1)) then
                print*,"cand: 1x",xe(i_range(i))
                dist_x = abs(xrange(in2) - xe(i_range(i)))
                dist_x2 = xrange(in1) - xe(i_range(i))
                dist_y = dist_x * tan(angle)
                if (abs(dist_x2).gt.tol .and. abs(dist_x).gt.tol) then
                  call AddToList_verts(intersect_p,(/xe(i_range(i)),yrange(1)+dist_y/))
                  print*, "what was added: ",(/xe(i_range(i)),yrange(1)+dist_y/)
                  no_int_par(k) = no_int_par(k) + 1
                  ci = ci + 1
                end if
              end if
            end do

            do i=1,size(j_range,1)
              print*,"cand:1 y",ye(j_range(i))

              if (ye(j_range(i)).lt.yrange(2) .and. ye(j_range(i)).gt.yrange(1)) then
                dist_y = abs(yrange(in2) - ye(j_range(i)))
                dist_y2 = yrange(in1) - ye(j_range(i))
                dist_x = dist_y/tan(angle)
                if (abs(dist_y2).gt.tol .and. abs(dist_y).gt.tol) then
                  call AddToList_verts(intersect_p,(/xrange(1)+dist_x,ye(j_range(i))/))
                  print*, "waht was added: ",(/xrange(1)+dist_x,ye(j_range(i))/)
                  no_int_par(k) = no_int_par(k) + 1
                  ci = ci + 1
                end if
              end if
            end do

            if (ci .gt. 1) then
              do i = 1,ci
                call AddToList_verts(pair_inter_par1,intersect_p(:,size(intersect_p,2)-(i-1)))
              end do
            end if
            if ( slope_bar.lt.0.d0) then
              i_0_p = i_0_p + logical_to_int(xe(i_0_p+1).lt.ver_hbox(1,2*k+1))
              j_0_p = j_0_p + sign(1.d0,slope_bar) * logical_to_int(ye(j_0_p).gt.ver_hbox(2,2*k+1))
            else if (slope_bar.gt.0.d0) then
              i_0_p = i_0_p + logical_to_int(xe(i_0_p+1).lt.ver_hbox(1,2*k+1)) &
                 - logical_to_int(xe(i_0_p).gt.ver_hbox(1,2*k+1))
              j_0_p = j_0_p + sign(1.d0,slope_bar)*logical_to_int(ye(j_0_p+1).lt.ver_hbox(2,2*k+1))&
              - sign(1.d0,slope_bar)*logical_to_int(ye(j_0_p).gt.ver_hbox(2,2*k+1))
            end if

            i_range2 = i_1 + i_chk_2
            j_range2 = j_1 + j_chk2

            ci=0
            do i = 1,size(i_range2,1)
              print*,"cand:2 x", xe(i_range2(i))
              print *, "Xrange: ",xrange2
! wrong here
              if (xe(i_range2(i)).lt.xrange2(2) .and. xe(i_range2(i)).gt.xrange2(1)) then
                dist_x2 = abs(xrange2(in2) - xe(i_range2(i)))
                dist_y2 = dist_x2 * tan(angle)
                dist_x3 = xrange2(in1) - xe(i_range2(i))
                if (abs(dist_x3).gt.tol .and. abs(dist_x2).gt.tol) then
                  call AddToList_verts(intersect_p,(/xe(i_range2(i)),yrange2(1)+dist_y2/))
                  print*, "waht was added: ",(/xe(i_range2(i)),yrange2(1)+dist_y2/)
                  no_int_par(k) = no_int_par(k) + 1
                  ci = ci+1
                end if
              end if
            enddo

            do i =1,size(j_range2,1)
              print*,"cand: 2y", ye(j_range2(i))

              if (ye(j_range2(i)).lt.yrange2(2) .and. ye(j_range2(i)).gt.yrange2(1)) then
                dist_y2 = abs(yrange2(in2) - ye(j_range2(i)))
                dist_x2 = dist_y2 / tan(angle)
                dist_y3 = yrange2(in1) - ye(j_range2(i))
                if (abs(dist_y3).gt.tol .and. abs(dist_y2).gt.tol) then
                  call AddToList_verts(intersect_p,(/xrange2(1)+dist_x2,ye(j_range2(i))/))
                  print*, "waht was added: ",(/xrange2(1)+dist_x2,ye(j_range2(i))/)
                  no_int_par(k) = no_int_par(k) + 1
                  ci= ci+1
                end if
              end if
            end do

            if (ci .gt. 1) then
              do i = 1,ci
                call AddToList_verts(pair_inter_par2,intersect_p(:,size(intersect_p,2)-(i-1)))
              end do
            end if

            if (k.eq.m)then
              exit
            end if
            if (slope_bar.lt.0.d0) then
              i_1 = i_1 + logical_to_int(xe(i_1+1).lt.ver_hbox(1,2*k+2))
              j_1 = j_1 + sign(1.d0,slope_bar)*logical_to_int(ye(j_1).gt.ver_hbox(2,2*k+2))
            else if (slope_bar.gt.0.d0) then
              i_1 =  i_1 + logical_to_int(xe(i_1+1).lt.ver_hbox(1,2*k+2)) &
                  - logical_to_int(xe(i_1).gt.ver_hbox(1,2*k+2))
              j_1 = j_1 + sign(1.d0,slope_bar)*logical_to_int(ye(j_1+1).lt.ver_hbox(2,2*k+2))&
               - sign(1.d0,slope_bar)*logical_to_int(ye(j_1).gt.ver_hbox(2,2*k+2))
             end if
           end do

          ! physical vertices inside or not:
          if (slope_bar.lt.0.d0) then
            i_h = i_0 + logical_to_int(xe(i_0+1).lt.hbox_cen(1,1))
            j_h = j_0 - logical_to_int(ye(j_0).gt.hbox_cen(2,1)) + &
              logical_to_int(ye(j_0+1).lt.hbox_cen(2,1))
          end if
          if (slope_bar.gt.0.d0) THEN
            i_h = i_0 + logical_to_int(xe(i_0+1).lt.hbox_cen(1,1)) &
            -logical_to_int(xe(i_0).gt.hbox_cen(1,1))
            j_h = j_0 - logical_to_int(ye(j_0).gt.hbox_cen(2,1)) + &
              logical_to_int(ye(j_0+1).lt.hbox_cen(2,1))
          end if
          call IAddToList_verts(rough_ind,(/i_h,j_h/))

          allocate(phys_ver_in(m))
          print *, "m:",m
          print*, "pair inter p:", pair_inter_par1
          print*, "Ih: ",i_h
          print*, "Xe(ih): ", xe(i_h)

          do k =1,m
            cands(1,:) = (/xe(i_h),xe(i_h),xe(i_h+1),xe(i_h+1)/)
            cands(2,:) = (/ye(j_h),ye(j_h+1),ye(j_h),ye(j_h+1)/)
            cen_x = hbox_cen(1,k)
            cen_y = hbox_cen(2,k)
            !hbox corners
            A = ver_hbox(:,2*k-1)
            B = ver_hbox(:,2*k)
            C = ver_hbox(:,2*k+1)
            D = ver_hbox(:,2*k+2)

            chk_in=0
            sides(1,:,:) = reshape((/A,B/),(/2,2/))
            sides(2,:,:) = reshape((/A,C/),(/2,2/))
            sides(3,:,:) = reshape((/C,D/),(/2,2/))
            sides(4,:,:) = reshape((/B,D/),(/2,2/))

            do i=1,4
              P = (/cands(1,i),cands(2,i)/)
              code = in_or_out(P,sides,4)
              if (code.eq.1 .and. dist(A,P).gt.tol .and. dist(B,P).gt.tol &
              .and. dist(C,P).gt.tol .and. dist(D,P).gt.tol) then
                chk_in = 1
                call AddToList_verts(phys_ver,P)
              end if
            end do
            if (chk_in.eq.1) then
              phys_ver_in(k) = 1
            else
              phys_ver_in(k) = 0
            end if
            if (k.eq.m) then
              exit
            end if
            if (slope_bar.lt.0.d0) then
              i_h = i_h + logical_to_int(xe(i_h+1).lt.hbox_cen(1,k+1))
              j_h = j_h - logical_to_int(ye(j_h).gt.hbox_cen(2,k+1))
            else
              i_h = i_h - logical_to_int(xe(i_h).gt.hbox_cen(1,k+1)) + &
                logical_to_int(xe(i_h+1).lt.hbox_cen(1,k+1))
              j_h = j_h + logical_to_int(ye(j_h+1).lt.hbox_cen(2,k+1)) -&
                logical_to_int(ye(j_h).gt.hbox_cen(2,k+1))
            end if
            call IAddToList_verts(rough_ind,(/i_h,j_h/))
          end do
          print*, "phys ver:", phys_ver

          ! Getting edges of hboxes
          if (allocated(pair_inter_norm)) then
            pair_inter_n = remove_dups(pair_inter_norm)
          endif
          if (allocated(pair_inter_par1)) then
            pair_inter_p1 = remove_dups(pair_inter_par1)
          end if
          if (allocated(pair_inter_par2)) then
            pair_inter_p2 = remove_dups(pair_inter_par2)
          end if
          allocate(no_int_norm2(m))
          no_int_norm2 = no_int_norm(2:) + no_int_norm(:m)

          print * ,"num norm",no_int_norm2
          print *, "num par:",no_int_par
          ! index counter/pointers
          i = 0
          j = 0
          j2 = 0
          k1 = 0
          k2 = 0
          k3 = 0
          ! num of edges in hbox
          allocate(no_edges(m))
          ! addending edges
          do k = 1,m  !!!!!!!!!!!!!!!!!!!!!!!!!! ID 0
            num_edges = 0
            num_norm = no_int_norm2(k)
            num_par = no_int_par(k)
            print*, "NUM Par: ", num_par
            A = ver_hbox(:,2*k-1)
            B = ver_hbox(:,2*k)
            C = ver_hbox(:,2*k+1)
            D = ver_hbox(:,2*k+2)
            x_vers = (/A(1),B(1),C(1),D(1)/)
            y_vers = (/A(2),B(2),C(2),D(2)/)
            call KB07AD(x_vers,4,dummy)
            call KB07AD(y_vers,4,dummy)
            xmin = x_vers(1)
            xmax = x_vers(4)
            ymin = y_vers(1)
            ymax = y_vers(4)

            normal_slope = ((C(2)-D(2))/(C(1)-D(1)))
            c_int = C(2) - normal_slope*C(1)
            do i2 = j,max(j+num_norm-1,0)
              if (i2.gt.size(intersect_n,2).or. num_norm.eq.0) then
                exit
              end if
              call AddToList_verts(inter_chk,intersect_n(:,i2+1))
              print*,"the intersection is: n ",intersect_n(:,i2+1)
              call AddToList_verts(inter_chk_n,intersect_n(:,i2+1))
            end do

            do i2 = j2, max(j2+num_par-1,0)
              if (i2 .gt. size(intersect_p,2) .or. num_par.eq.0) then
                exit
              end if
              call AddToList_verts(inter_chk,intersect_p(:,i2+1))
              print*,"the intersection is: p ", intersect_p(:,i2+1)
            end do
            if (.not. allocated(inter_chk_n)) then
              j = j + num_norm
            else if (allocated(inter_chk_n)) then
              if (size(inter_chk_n,2).eq.1) then
                last_ver = inter_chk_n(:,1)
                if (abs(last_ver(2)-(normal_slope*last_ver(1) + c_int)).lt.tol) then
                  j = j + num_norm - 1
                else
                  j = j + 1
                end if
              end if
              if (size(inter_chk_n,2).gt.1) then
                last_ver = inter_chk_n(:,size(inter_chk_n,2))
                pen_ver = inter_chk_n(:,size(inter_chk_n,2)-1)
                if (abs(last_ver(2)-(normal_slope*last_ver(1)+c_int)).lt.tol &
                .and. abs(pen_ver(2)-(normal_slope*pen_ver(1) +c_int)).lt.tol) then
                  j = j + num_norm - 2
                else if (abs(last_ver(2)-(normal_slope*last_ver(1)+c_int)).lt.tol) then
                  j = j +num_norm -1
                else
                  j = j +num_norm
                end if
              end if
            end if

            j2 = j2 + num_par

            ! remove duplicates:
            if (allocated(inter_chk)) then
              inter_chkd = remove_dups(inter_chk)
            end if
            print*,"intersections:", inter_chk

            ! get all the right edges:
            if (.not. allocated(inter_chk)) then
              call AddToList_edges(edges,reshape((/A,B/),(/2,2/)))
              call AddToList_edges(edges,reshape((/A,C/),(/2,2/)))
              call AddToList_edges(edges,reshape((/C,D/),(/2,2/)))
              call AddToList_edges(edges,reshape((/B,D/),(/2,2/)))
              num_edges = num_edges + 4
            else
              allocate(wherea(size(inter_chkd,2),4))
              wherea = 0
              if (phys_ver_in(k).eq.1) then
                if (i.ge.size(phys_ver,2)) then
                  i = i - 1
                end if
                x_phys_ver = phys_ver(1,i+1)
                y_phys_ver = phys_ver(2,i+1)
                print *, "x,y", x_phys_ver,y_phys_ver

                i = i + 1
                if (i.ge.size(phys_ver,2)) then
                  i = i - 1
                end if
                call AddToList_verts(phizverts,(/x_phys_ver,y_phys_ver/))
                sides(1,:,:) = reshape((/A,B/),(/2,2/))
                sides(2,:,:) = reshape((/A,C/),(/2,2/))
                sides(3,:,:) = reshape((/C,D/),(/2,2/))
                sides(4,:,:) = reshape((/B,D/),(/2,2/))
                P = (/phys_ver(1,i+1),phys_ver(2,i+1)/)
                code = in_or_out(P,sides,4)
                if (code.eq.1) then
                  call AddToList_verts(phizverts,(/phys_ver(1,i+1),phys_ver(2,i+1)/))
                  i = i+1
                end if
                  d_check = huge(0.d0)
                if (allocated(phizverts)) then         !! need an end if statement somewhere down below
                  if (size(phizverts,2).gt.1) then
                    phiz2 = phizverts(:,2)
                    phiz1 = phizverts(:,1)
                    d_check = dist(phiz1,phiz2)
                    if (d_check.lt.tol) then
                      d_check = huge(0.d0)
                      call DelFromList_verts(phizverts)
                    end if
                   end if
                   do p1=1,size(phizverts,2)  !! HERE 1
                     x_phiz = phizverts(1,p1)
                     y_phiz = phizverts(2,p1)
                     do m1 = 1,size(inter_chkd,2) !! HERE 2
                        if (inter_chkd(1,m1) .gt. xmax .or. inter_chkd(1,m1).lt.xmin .or. &
                        inter_chkd(2,m1).gt. ymax .or. inter_chkd(2,m1).lt.ymin) then
                        continue
                  elseif ((abs(x_phiz-inter_chkd(1,m1)).lt.tol.or.abs(y_phiz-inter_chkd(2,m1)).lt.tol)&
                        .and.(abs(x_phiz-inter_chkd(1,m1)).gt.tol.or.abs(y_phiz-inter_chkd(2,m1)).gt.tol)) then
                    if (size(phizverts,2).gt.1) then
                      if (dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).lt.d_check .and. &
                      dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).gt.tol) then
                        v1 = (/x_phiz,y_phiz/)
                        call AddToList_edges(edges,reshape((/v1,inter_chkd(:,m1)/),(/2,2/)))
                        num_edges = num_edges + 1
                      end if
                    else if (dist((/x_phiz,y_phiz/),inter_chkd(:,m1)).gt.tol) then
                      v1 = (/x_phiz,y_phiz/)
                      call AddToList_edges(edges,reshape((/v1,inter_chkd(:,m1)/),(/2,2/)))
                      num_edges = num_edges + 1
                    endif
                  end if
                  normal_slope1 = (A(2)-B(2))/(A(1)-B(1))
                  c1 = A(2) - normal_slope1*A(1)
                  normal_slope2 = (A(2)-C(2))/(A(1)-C(1))
                  c2 = A(2) - normal_slope2*A(1)
                  normal_slope3 = (D(2)-B(2))/(D(1)-B(1))
                  c3 = D(2) - normal_slope3*D(1)
                  normal_slope4 = (D(2)-C(2))/(D(1)-C(1))
                  c4 = D(2) - normal_slope4*D(1)
                  if (abs(inter_chkd(2,m1)-(normal_slope1*inter_chkd(1,m1)+c1)).lt.tol) then
                    wherea(m1,1) = wherea(m1,1) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope2*inter_chkd(1,m1)+c2)).lt.tol) then
                    wherea(m1,2) = wherea(m1,2) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope3*inter_chkd(1,m1)+c3)).lt.tol) then
                    wherea(m1,3) = wherea(m1,3) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope4*inter_chkd(1,m1)+c4)).lt.tol) then
                    wherea(m1,4) = wherea(m1,4) + 1
                  end if
                end do

                ! get which row has nonzero for each column, thats the index of the intersections.
                print *, "PHIZ:",x_phiz,y_phiz
                print *, " A:", A
                print *, " B: ",B
                print *, "C:", C
                print *, " D:",D
                if ((abs(x_phiz-A(1)).lt.tol .or. abs(y_phiz-A(2)).lt.tol) .and. &
                   dist(phizverts(:,p1),A).lt.d_check ) then
                   v1 = (/x_phiz,y_phiz/)
                   call AddToList_edges(edges,reshape((/v1,A/),(/2,2/)))
                   num_edges = num_edges + 1
                 end if

                 if ((abs(x_phiz-B(1)).lt.tol .or. abs(y_phiz-B(2)).lt.tol) .and. &
                    dist(phizverts(:,p1),B).lt.d_check ) then
                    v1 = (/x_phiz,y_phiz/)
                    call AddToList_edges(edges,reshape((/v1,B/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if

                 if ((abs(x_phiz-C(1)).lt.tol .or. abs(y_phiz-C(2)).lt.tol) .and. &
                    dist(phizverts(:,p1),C).lt.d_check ) then
                    v1 = (/x_phiz,y_phiz/)
                    call AddToList_edges(edges,reshape((/v1,C/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if

                 if ((abs(x_phiz-D(1)).lt.tol .or. abs(y_phiz-D(2)).lt.tol) .and. &
                    dist(phizverts(:,p1),D).lt.d_check ) then
                    v1 = (/x_phiz,y_phiz/)
                    call AddToList_edges(edges,reshape((/v1,D/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                end do
              end if
              print*, "phiz vert", phizverts
              ! print *, "interchkd", inter_chkd
              print *, "size intchkd", size(inter_chkd,2)
                if (allocated(phizverts)) then  !#################   ID 1 need end if somewhere down
                  if (size(phizverts,2).eq.1) then   !! ID 2
                    phiz = phizverts(:,1)
                    x_phiz = phiz(1)
                    y_phiz = phiz(2)

                    allocate(listt(size(inter_chkd,2)))
                    do m1 = 1,size(inter_chkd,2)
                      listt = (/(i2,i2=1,size(inter_chkd,2))/)
                      call IDelFromList(listt,m1)
                      do m2 = 1,size(listt,1)
                        m3 = listt(m2)
                        if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                        abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                        (abs(inter_chkd(1,m1)-x_phiz).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz).gt.tol)&
                        .and.(abs(inter_chkd(1,m1)-inter_chkd(1,m3)).gt.tol.or.abs(inter_chkd(2,m1)-&
                        inter_chkd(2,m3)).gt.tol)) then
                          call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                          num_edges = num_edges + 1
                        end if
                      end do
                    end do
                  ! else if (size(phizverts,2).gt.1) then   !! ID 2 IF CLAUSE
                  !   phiz1 = phizverts(:,1)
                  !   x_phiz1 = phiz1(1)
                  !   y_phiz1 = phiz1(2)
                  !   phiz2 = phizverts(:,2)
                  !   x_phiz2 = phiz2(1)
                  !   y_phiz2 = phiz2(2)
                  !   do m1 = 1, size(inter_chkd,2)
                  !     listt = (/(i2,i2=1,size(inter_chkd,2))/)
                  !     call IDelFromList(listt,m1)
                  !     do m2 = 1,size(listt,1)
                  !       m3 = listt(m2)
                  !       if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                  !       abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                  !       (abs(inter_chkd(1,m1)-x_phiz1).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz1).gt.tol)&
                  !  .and.(abs(inter_chkd(1,m1)-x_phiz2).gt.tol.and.abs(inter_chkd(2,m1)-y_phiz2).gt.tol)&
                  !       .and.(abs(inter_chkd(1,m1)-inter_chkd(1,m3)).gt.tol.or.abs(inter_chkd(2,m1)-&
                  !       inter_chkd(2,m3)).gt.tol)) then
                  !       call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                  !       num_edges = num_edges + 1
                  !       end if
                  !     end do
                  !   end do
                  !   if ((abs(x_phiz1-x_phiz2).lt.tol.or.abs(y_phiz1-y_phiz2).lt.tol).and. &
                  !      dist(phiz1,phiz2).gt.tol) then
                  !      m3 = size(phizverts,2)
                  !      call AddToList_edges(edges,reshape((/phizverts(:,1),phizverts(:,m3)/),(/2,2/)))
                  !      num_edges = num_edges + 1
                  !   end if
                  end if
                end if

              else  !!!!  ID 3
                do m1 = 1,size(inter_chkd,2)
                  listt = (/(i2,i2=1,size(inter_chkd,2))/)
                  call IDelFromList(listt,m1)
                  do m2 = 1,size(listt,1)
                    m3 = listt(m2)
                    if ((abs(inter_chkd(1,m1)-inter_chkd(1,m3)).lt.tol .or. &
                    abs(inter_chkd(2,m1)-inter_chkd(2,m3)).lt.tol) .and. &
                       dist(inter_chkd(:,m1),inter_chkd(:,m3)).gt.tol) then
                      call AddToList_edges(edges,reshape((/inter_chkd(:,m1),inter_chkd(:,m3)/),(/2,2/)))
                      num_edges = num_edges + 1
                    end if
                  end do
                  normal_slope1 = (A(2)-B(2))/(A(1)-B(1))
                  c1 = A(2) - normal_slope1*A(1)
                  normal_slope2 = (A(2)-C(2))/(A(1)-C(1))
                  c2 = A(2) - normal_slope2*A(1)
                  normal_slope3 = (D(2)-B(2))/(D(1)-B(1))
                  c3 = D(2) - normal_slope3*D(1)
                  normal_slope4 = (D(2)-C(2))/(D(1)-C(1))
                  c4 = D(2) - normal_slope4*D(1)
                  if (abs(inter_chkd(2,m1)-(normal_slope1*inter_chkd(1,m1)+c1)).lt.tol) then
                    wherea(m1,1) = wherea(m1,1) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope2*inter_chkd(1,m1)+c2)).lt.tol) then
                    wherea(m1,2) = wherea(m1,2) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope3*inter_chkd(1,m1)+c3)).lt.tol) then
                    wherea(m1,3) = wherea(m1,3) + 1
                  else if (abs(inter_chkd(2,m1)-(normal_slope4*inter_chkd(1,m1)+c4)).lt.tol) then
                    wherea(m1,4) = wherea(m1,4) + 1
                  end if
                end do
              end if

              x_cat(1,:) = (/A(1),A(1),B(1),C(1)/)
              x_cat(2,:) = (/B(1),C(1),D(1),D(1)/)
              y_cat(1,:) = (/A(2),A(2),B(2),C(2)/)
              y_cat(2,:) = (/B(2),C(2),D(2),D(2)/)
              where_sum = sum(wherea,dim=1)

              do m3 = 1,4
                v1 = (/x_cat(1,m3),y_cat(1,m3)/)
                v2 = (/x_cat(2,m3),y_cat(2,m3)/)
                if ( where_sum(m3) .eq. 0 .and.dist(v1,v2).gt.tol) then
                 call AddToList_edges(edges,reshape((/v1,v2/),(/2,2/)))
                 num_edges = num_edges + 1
                end if
                if (where_sum(m3) .eq. 1) then
                  o = nonzero(wherea(:,m3))
                  if (dist(inter_chkd(:,o),v1).gt.tol)THEN
                    call AddToList_edges(edges,reshape((/inter_chkd(:,o),v1/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                  if (dist(inter_chkd(:,o),v2).gt.tol)then
                    call AddToList_edges(edges,reshape((/inter_chkd(:,o),v2/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                end if
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1333
                if (where_sum(m3) .gt. 1) then
                  find = nonzero(wherea(:,m3))
                  allocate(x_list(size(find)))
                  x_list = inter_chkd(1,find)
                  min_loc = minloc(x_list)
                  ml = min_loc(1)
                  max_loc = maxloc(x_list)
                  nl = max_loc(1)
                  ver_x = x_cat(:,m3)
                  ver_min = minloc(ver_x)
                  vl = ver_min(1)
                  ver_max = maxloc(ver_x)
                  wl = ver_max(1)
                  v3 = inter_chkd(:,find(ml))
                  v4 = (/minval(x_cat(:,m3)),y_cat(vl,m3)/)
                  if (dist(v3,v4).gt.tol) then
                    call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                  v3 = inter_chkd(:,find(nl))
                  v4 = (/maxval(x_cat(:,m3)),y_cat(wl,m3)/)
                  if (dist(v3,v4).gt.tol) then
                    call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                  v4 = inter_chkd(:,find(ml))
                  if (dist(v3,v4).gt.tol) then
                    call AddToList_edges(edges,reshape((/v3,v4/),(/2,2/)))
                    num_edges = num_edges + 1
                  end if
                  deallocate(x_list)
                end if
              end do
            end if

            no_edges(k) =  num_edges
            if (allocated(wherea)) then
              deallocate(wherea)
            end if
            if (allocated(phizverts)) then
              deallocate(phizverts)
            end if
            if (allocated(listt)) then
              deallocate(listt)
            end if
            if (allocated(inter_chkd)) then
              deallocate(inter_chkd)
              deallocate(inter_chk)
            end if
          end do

          !!!!!!!!!!!1! getting the fragments of each hbox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          counter = 0
          counter2 = 0
          ! allocate(q_hbox(3,m))
          ! allocate(aux_hbox(1,m))
          ! q_hbox = 0.d0
          ! aux_hbox = 0.d0
          do k = 1, m
            print *,"K:",k
            num_edge = no_edges(k)
            print*, "num edges",num_edge
            do j=counter,max(counter+num_edge-1,0)
              e1 = edges(:,:,j+1)
              ord = vec2ord(e1(1,:))
              call AddToList_edges(edges_chk,e1(:,ord))
            end do
            ! deallocate(indexes)
            ! allocate(indexes(size(edges_chk,3)))
            ! call KB07AD(edges_chk(1,1,:),size(edges_chk,3),indexes)
            ! edges_chk(:,:,1:size(edges_chk,3)) = edges_chk(:,:,indexes)
            counter = counter + num_edge
            print* , "edges:",edges_chk
            ! hbox area info and list to collect areas and indices of frags
            ! hbox vertices
            A = ver_hbox(:,2*k-1)
            B = ver_hbox(:,2*k)
            C = ver_hbox(:,2*k+1)
            D = ver_hbox(:,2*k+2)
            hbox_area = area_polygon((/A(1),B(1),D(1),C(1),A(1)/),(/A(2),B(2),D(2),C(2),A(2)/) )!dist(A,B)*dist(B,D)
            print *, "HBOX AREA: ",hbox_area
            init_area = 0.d0

            if (num_edge .eq. 4) then
              call AddToList_verts(cyclea,A)
              call AddToList_verts(cyclea,B)
              call AddToList_verts(cyclea,D)
              call AddToList_verts(cyclea,C)
              call AddToList_verts(cyclea,A)
              call AddToList(areas,hbox_area)
              call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
              ,xe,ye,mx,my,mbc,right_ind,centroid)
              call IAddToList_verts(inds,right_ind)
              init_area = init_area + hbox_area
              num_frags(k) = 1
              hbox_areas(k) = hbox_area
              area_frags(1,k) = 1.d0
              index_frags(:,1,k) = right_ind
              goto 99
            end if
            if (allocated(cyclea)) then
              deallocate(cyclea)
            end if
            ! get the lateral two edges, "start edges" to do the rotation
            do j = 1, size(edges_chk,3)
              e1 = edges_chk(:,:,j)
              v4 = e1(:,2)
              v5 = e1(:,1)
              v6 = v4 - v5
              if (abs(v6(2)) .lt. tol) then
                if (.not. allocated(horz_edges)) then
                  call AddToList_edges(horz_edges,e1)
                else ! check if edge e1 is not in horz edges and e1 permuted not in horz edges
                  if (.not. in(e1,horz_edges)) then
                    call AddToList_edges(horz_edges,e1)
                  end if
                end if
              end if
              if (abs(v6(1)).lt.tol) then
                if (.not. allocated(vert_edges)) then
                  call AddToList_edges(vert_edges,e1)
                else
                  if (.not. in(e1,vert_edges)) then
                    call AddToList_edges(vert_edges,e1)
                  end if
                end if
              end if
            end do
            ! if both horz and vert edges exist
            print*, "vert edge ", vert_edges
            if (allocated(horz_edges)) then
              print*, "horze edge", horz_edges
              ! horz edge counter/clkwise cycle search
              isign = -1
              do j=1,size(horz_edges,3)
                h_edge = horz_edges(:,:,j)
                ! print*, "here"
                do i = 1,2
                  ! print *, "region", region_d
                  ! print *, "cycle", cyclea
                  call AddToList_edges(region_d,h_edge)
                  call AddToList_verts(cyclea,h_edge(:,1))
                  call AddToList_verts(cyclea,h_edge(:,2))
                  !!!!!!!!!!! here 9.18.20 5:38 PM
                  ! print*, "stop"

                  call is_cycle(reshape(cyclea,(/2,2,1/)),yes_or_no)
                  i2 = 1

                  do while (.not.yes_or_no)
                    do i3=1,size(edges_chk,3)
                      e1 = edges_chk(:,:,i3)
                      v1 = e1(:,1)
                      v2 = e1(:,2)
                      x1 = v1(1)
                      x2 = v2(1)
                      y1 = v1(2)
                      y2 = v2(2)
                      last_elt = region_d(:,:,size(region_d,3))
                      v_last1 = last_elt(:,1)
                      v_last2 = last_elt(:,2)

                      if ((dist(v1,v_last2).lt.tol .or. dist(v2,v_last2).lt.tol).or.&
                      (dist(v1,v_last1).lt.tol .or. dist(v2,v_last1).lt.tol)) then
                         if (dist(v1,v_last2).lt.tol) then
                            cand = v2
                          else if (dist(v2,v_last2).lt.tol) then
                            cand = v1
                          else if (dist(v1,v_last1).lt.tol) then
                            cand = v2
                          else if (dist(v2,v_last1).lt.tol) then
                            cand = v1
                          end if

                          ver3 = cand
                          ver2 = cyclea(:,i2+1)
                          ver1 = cyclea(:,i2)
                          vec1 = ver3 - ver2
                          vec2 = ver2 - ver1

                          z = vec2(1) * vec1(2) - vec1(1)*vec2(2)
                          if (sign(1.d0,z*isign).gt.0.d0 .and. abs(z) .gt.tol**2) then
                            if (.not. in(e1, region_d).and..not.vert_in_edge(cyclea(:,i2),e1))&
                            then
                              call AddToList_verts(cyclea,cand)
                              call AddToList_edges(region_d,e1)
                              call is_cycle(region_d,yes_or_no)
                              i2 = i2 + 1
                              exit
                            end if
                          endif
                          continue
                        end if
                      end do
                      call is_cycle(region_d,yes_or_no)
                    end do ! while
                    area = area_polygon(cyclea(1,:),cyclea(2,:))
                    print*, "Area ", area
                    print*, "cycle:", cyclea
                    call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
                     ,xe,ye,mx,my,mbc,right_ind,centroid)
                    print*, " right ind: ", right_ind
                    if (.not. allocated(areas)) then
                      init_area = init_area + abs(area)
                      call AddToList(areas,area)
                      call IAddToList_verts(inds,right_ind)
                    else if (.not. ind_in_inds(right_ind,inds)) THEN
                      init_area= init_area + area
                      call AddToList(areas,area)
                      call IAddToList_verts(inds,right_ind)
                      if (abs(init_area-hbox_area).lt.100.d0*tol) then
                        goto 99
                      end if
                    end if

                    isign = isign * -1
                    ! print*, "here 3"

                    if (allocated(region_d)) then
                      deallocate(region_d)
                    end if
                    if (allocated(cyclea)) then
                      deallocate(cyclea)
                    end if

                  end do
                end do
              end if

              ! VERT EDGES :
              if (allocated(vert_edges)) then
                print*, "vert edge ", vert_edges
                ! horz edge counter/clkwise cycle search
                isign = -1
                do j=1,size(vert_edges,3)
                  v_edge = vert_edges(:,:,j)
                  do i = 1,2
                    call AddToList_edges(region_d,v_edge)
                    call AddToList_verts(cyclea,v_edge(:,1))
                    call AddToList_verts(cyclea,v_edge(:,2))
                    !!!!!!!!!!! here 9.18.20 5:38 PM
                    call is_cycle(reshape(cyclea,(/2,2,1/)),yes_or_no)
                    i2 = 1
                    do while (.not.yes_or_no)
                      do i3=1,size(edges_chk,3)
                        e1 = edges_chk(:,:,i3)
                        v1 = e1(:,1)
                        v2 = e1(:,2)
                        x1 = v1(1)
                        x2 = v2(1)
                        y1 = v1(2)
                        y2 = v2(2)
                        last_elt = region_d(:,:,size(region_d,3))
                        v_last1 = last_elt(:,1)
                        v_last2 = last_elt(:,2)

                        if ((dist(v1,v_last2).lt.tol .or. dist(v2,v_last2).lt.tol).or.&
                        (dist(v1,v_last1).lt.tol .or. dist(v2,v_last1).lt.tol)) then
                           if (dist(v1,v_last2).lt.tol) then
                              cand = v2
                            else if (dist(v2,v_last2).lt.tol) then
                              cand = v1
                            else if (dist(v1,v_last1).lt.tol) then
                              cand = v2
                            else if (dist(v2,v_last1).lt.tol) then
                              cand = v1
                            end if

                            ver3 = cand
                            ver2 = cyclea(:,i2+1)
                            ver1 = cyclea(:,i2)
                            vec1 = ver3 - ver2
                            vec2 = ver2 - ver1

                            z = vec2(1) * vec1(2) - vec1(1)*vec2(2)
                            if (sign(1.d0,z*isign).gt.0.d0 .and. abs(z) .gt.tol**2 ) then
                              if (.not. in(e1, region_d) .and. .not. vert_in_edge(cyclea(:,i2),e1))&
                              then
                                call AddToList_verts(cyclea,cand)
                                call AddToList_edges(region_d,e1)
                                call is_cycle(region_d,yes_or_no)
                                i2 = i2 + 1
                                exit
                              end if
                            endif
                            continue
                          end if
                        end do
                      end do ! while
                      area = area_polygon(cyclea(1,:),cyclea(2,:))
                      print *, " AREA" , area
                      print*, cyclea
                      call find_ind(cyclea,rough_ind(:,min(k,size(rough_ind,2)))&
                       ,xe,ye,mx,my,mbc,right_ind,centroid)
                       print *, "RI",right_ind
                      if (.not. allocated(areas)) then
                        init_area = init_area + area
                        call AddToList(areas,area)
                        call IAddToList_verts(inds,right_ind)
                      else if (.not. ind_in_inds(right_ind,inds)) THEN
                        init_area= init_area + area
                        call AddToList(areas,area)
                        call IAddToList_verts(inds,right_ind)
                        if (abs(init_area-hbox_area).lt.tol*100.d0) then
                          goto 99
                        end if
                      end if
                      isign = isign * -1
                      deallocate(region_d)
                      deallocate(cyclea)
                    end do
                  end do
                end if
      99 continue
                print *, "Areas: ", init_area, hbox_area
                num_frags(k) = size(areas,1)
                hbox_areas(k) = hbox_area
                if (abs(init_area-hbox_area).gt.tol) THEN
                  print *, "your code didn't work!!"
                end if
                print*, "INDS",inds
                do i=1,size(areas,1)
                   look = inds(:,i)
                   area_frags(i,k) = areas(i)/hbox_area
                   index_frags(:,i,k) = look
                   ! q_hbox(:,k) = q_hbox(:,k) + q_old(:,look(1),look(2))*areas(i)/hbox_area
                   ! aux_hbox(k) = aux_hbox(k) + aux(1,look(1),look(2))*areas(i)/hbox_area
                 enddo
                 deallocate(edges_chk)
                 if (allocated(horz_edges)) then
                   deallocate(horz_edges)
                 end if
                 if (allocated(vert_edges)) then
                   deallocate(vert_edges)
                 end if
                 if (allocated(areas)) then
                   deallocate(areas)
                 end if
                 if (allocated(inds)) then
                   deallocate(inds)
                 end if
                 if (allocated(region_d)) then
                   deallocate(region_d)
                 end if
                 if (allocated(cyclea)) then
                   deallocate(cyclea)
                 end if
               end do

    end subroutine

    ! subroutine cut_cells_inds(hbox_cen_x,hbox_cen_y,mx,my,mbc,xe,ye,ii,jj)
    !   ! this subroutine finds the index-pairs of the cut cells given a straight line of embedded barrier
    !   !  input: x_0,y_0,x_e,y_e :: real(8) array of N intersecting points of barrier with grid
    !   !         i_0,j_0,i_e,j_e :: integer of start and end cells' indices of barrier
    !   !         x_lower,y_lower the lower left corner of grid
    !   !         dx,dy are mesh sizes of grid
    !   !         maxm is max(mx,my) the max number of grid cells in a row or col
    !   !         mbc is number of ghost cells
    !   ! output: ii the x-index of cut cells (1D array of integers)
    !   !         jj the y-index of cut cells (1D array of integers)
    !   !         num_cut_cells is the number of cut cells
    !
    !   implicit none
    !   real(kind=8) :: hbox_cen_x(:),hbox_cen_y(:)
    !   integer :: i_0,j_0,i,m
    !   integer :: mbc,mx,my
    !   integer,allocatable :: ii(:),jj(:)
    !   real(kind=8) :: xe(-mbc:mx+mbc), ye(-mbc:my+mbc)
    !
    !   m = len(hbox_cen_x)
    !   allocate(ii(m))
    !   allocate(jj(m))
    !   ! put them into ii and jj
    !   ! the first cell where the barrier starts is always cut:
    !   do i=-mbc,mx+mbc-1
    !     if (hbox_cen_x(1) .gt. xe(i) .and. hbox_cen_x(1) .lt. xe(i+1)) then
    !       i_0 = i+1
    !       exit
    !     end if
    !   end do
    !   do i = -mbc,my+mbc-1
    !     if (hbox_cen_y(1) .gt. ye(i) .and. hbox_cen_y(1) .lt. ye(i+1)) then
    !       j_0 = i+1
    !       exit
    !     end if
    !   end do
    !
    !   ii(1) = i_0
    !   jj(1) = j_0
    !   do i = 1, m-1
    !     if (hbox_cen_x(i+1) .lt. xe(ii(i)) .and. hbox_cen_x(i+1) .gt. xe(ii(i)-1)) then
    !       ii(i+1) = ii(i) - 1
    !     else if (hbox_cen_x(i+1) .gt. xe(ii(i)+1) .and. hbox_cen_x(i+1) .lt. xe(ii(i)+2)) then
    !       ii(i+1) = ii(i) + 1
    !     else if (hbox_cen_x(i+1) .gt. xe(ii(i)) .and. hbox_cen_x(i+1) .lt. xe(ii(i)+1)) then
    !       ii(i+1) = ii(i)
    !     end if
    !     if (hbox_cen_y(i+1) .lt. ye(jj(i)) .and. hbox_cen_y(i+1) .gt. ye(jj(i)-1)) then
    !       jj(i+1) = jj(i) - 1
    !     else if (hbox_cen_y(i+1) .gt. ye(jj(i)+1) .and. hbox_cen_y(i+1) .lt. ye(jj(i)+2)) then
    !       jj(i+1) = jj(i) + 1
    !     else if (hbox_cen_y(i+1) .gt. ye(jj(i)) .and. hbox_cen_y(i+1) .lt. ye(jj(i)+1)) then
    !       jj(i+1) = jj(i)
    !     end if
    !   end do
    !
    ! end subroutine cut_cells_inds








end module hbox_layer

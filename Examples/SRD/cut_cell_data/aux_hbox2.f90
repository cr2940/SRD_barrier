module aux_module_hbox

  use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
contains

  subroutine AddToList(list, element)

     IMPLICIT NONE

     integer :: i, isize
     double precision, intent(in) :: element
     double precision, dimension(:), allocatable, intent(inout) :: list
     double precision, dimension(:), allocatable :: clist


     if(allocated(list)) then
         isize = size(list)
         allocate(clist(isize+1))
         do i=1,isize
         clist(i) = list(i)
         end do
         clist(isize+1) = element

         deallocate(list)
         call move_alloc(clist, list)

     else
         allocate(list(1))
         list(1) = element
     end if


 end subroutine AddToList

      subroutine DelFromList_verts(list)
        ! delete last part of list
        IMPLICIT NONE

        integer :: i, isize
        double precision, dimension(:,:), allocatable, intent(inout) :: list
        double precision, dimension(:,:), allocatable :: clist


        if(allocated(list)) then
            isize = size(list,2)
            allocate(clist(2,isize-1))
            do i=1,isize-1
            clist(:,i) = list(:,i)
            end do

            deallocate(list)
            call move_alloc(clist, list)
        end if
      end subroutine

      subroutine DelFromList_edges(list,n)
        ! delete nth entry from list of edges
        implicit none
        integer :: i, isize, n , j
        real(8), allocatable , intent(inout):: list(:,:,:)
        real(8), allocatable :: clist(:,:,:)

        if (allocated(list)) then
          j = 1
          isize = size(list,3)
          allocate(clist(2,2,isize-1))
          do i = 1,isize
            if (i.ne.n) then
              clist(:,:,j) = list(:,:,i)
              j = j + 1
            end if
          end do
          deallocate(list)
          call move_alloc(clist, list)
        end if
      end subroutine

      subroutine DelFromList_nvert(list,n)
        implicit none
        integer :: i, isize, n , j
        real(8), allocatable , intent(inout):: list(:,:)
        real(8), allocatable :: clist(:,:)
        ! print*, "list: ",list
        ! print *, "n: ", n
        if (n .gt. size(list,2)) then
          return
        end if
        if (allocated(list)) then
          j = 1
          isize = size(list,2)
          allocate(clist(2,isize-1))
          do i = 1,isize
            if (i.ne.n) then
              clist(:,j) = list(:,i)
              j = j + 1
            end if
          end do
          deallocate(list)
          call move_alloc(clist, list)
        end if
      end subroutine



       ! adapted from https://stackoverflow.com/questions/28048508/how-to-add-new-element-to-dynamical-array-in-fortran90
       subroutine AddToList_verts(list, element)
       ! adding (/a,b/) to a list of two-elem vectors or empty list
         IMPLICIT NONE

         integer :: i, isize
         double precision, intent(in) :: element(2)
         double precision, dimension(:,:), allocatable, intent(inout) :: list
         double precision, dimension(:,:), allocatable :: clist


         if(allocated(list)) then
             isize = size(list,2)
             allocate(clist(2,isize+1))
             do i=1,isize
             clist(:,i) = list(:,i)
             end do
             clist(:,isize+1) = element

             deallocate(list)
             call move_alloc(clist, list)

         else
             allocate(list(2,1))
             list(:,1) = element
         end if
        end subroutine AddToList_verts


        ! adapted from https://stackoverflow.com/questions/28048508/how-to-add-new-element-to-dynamical-array-in-fortran90
        subroutine AddToList_edges(list, element)
        ! adding ((/a,b/),(/c,d/)) to a list of 2x2 matrices or empty list (represent edges)
          IMPLICIT NONE

          integer :: i, isize
          double precision, intent(in) :: element(2,2)
          double precision, dimension(:,:,:), allocatable, intent(inout) :: list
          double precision, dimension(:,:,:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list,3)
              allocate(clist(2,2,isize+1))
              do i=1,isize
              clist(:,:,i) = list(:,:,i)
              end do
              clist(:,:,isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(2,2,1))
              list(:,:,1) = element
          end if
        end subroutine AddToList_edges

        function nonzero(array) result(indexa)
          implicit none
          integer :: array(:), i,L
          integer, allocatable :: indexa(:)
          L = size(array,1)
          array = array.ne.0
          array = array*(/(i,i=1,L)/)
          indexa = pack(array,array/=0)
        end function

        function vec2ord(vec) result(ord)
          implicit none
          real(8) :: vec(2)
          integer :: ord(2)

          if (vec(1).le.vec(2)) then
            ord = (/1,2/)
          else
            ord = (/2,1/)
          end if
        end function

        function in(edge, edges_list) result(yes)
          implicit none
          real(8) :: edge(2,2)
          real(8),intent(in) :: edges_list(:,:,:)
          integer :: L , i
          logical :: yes
          yes = .false.
          L = size(edges_list,3)
          do i=1,L
            if (all(edge.eq.edges_list(:,:,i)) .or. &
              all(edge(:,(/2,1/)).eq.edges_list(:,:,i)))then
              yes = .true.
            end if
          end do
        end function

        subroutine IDelFromList(list, element)

           IMPLICIT NONE

           integer :: i, isize, j
           integer, intent(in) :: element
           integer, dimension(:), allocatable, intent(inout) :: list
           integer, dimension(:), allocatable :: clist

           j = 1
           if(allocated(list)) then
               isize = size(list)
               allocate(clist(isize-1))
               do i=1,isize
                 if (i .ne. element) then
                    clist(j) = list(i)
                    j = j + 1
                  end if
               end do

               deallocate(list)
               call move_alloc(clist, list)
           end if

        end subroutine IDelFromList






        subroutine IAddToList(list, element)

           IMPLICIT NONE

           integer :: i, isize
           integer, intent(in) :: element
           integer, dimension(:), allocatable, intent(inout) :: list
           integer, dimension(:), allocatable :: clist


           if(allocated(list)) then
               isize = size(list)
               allocate(clist(isize+1))
               do i=1,isize
               clist(i) = list(i)
               end do
               clist(isize+1) = element

               deallocate(list)
               call move_alloc(clist, list)

           else
               allocate(list(1))
               list(1) = element
           end if

       end subroutine IAddToList

        subroutine IAddToList_verts(list, element)
        ! adding (/a,b/) to a list of two-elem vectors or empty list
          IMPLICIT NONE

          integer :: i, isize
          integer, intent(in) :: element(2)
          integer, dimension(:,:), allocatable, intent(inout) :: list
          integer, dimension(:,:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list,2)
              allocate(clist(2,isize+1))
              do i=1,isize
              clist(:,i) = list(:,i)
              end do
              clist(:,isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(2,1))
              list(:,1) = element
          end if
        end subroutine IAddToList_verts

        subroutine IAddToList_edges(list, element)
        ! adding ((/a,b/),(/c,d/)) to a list of 2x2 matrices or empty list (represent edges)
          IMPLICIT NONE

          integer :: i, isize
          integer, intent(in) :: element(2,2)
          integer, dimension(:,:,:), allocatable, intent(inout) :: list
          integer, dimension(:,:,:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list,3)
              allocate(clist(2,2,isize+1))
              do i=1,isize
              clist(:,:,i) = list(:,:,i)
              end do
              clist(:,:,isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(2,2,1))
              list(:,:,1) = element
          end if
        end subroutine IAddToList_edges

        function dist(A,B) result(distance)
          implicit none
          real(kind=8) :: A(2),B(2)
          real(kind=8) :: distance
          distance = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
        end function dist

      ! below from : https://www.rosettacode.org/wiki/Remove_duplicate_elements#Fortran
        function remove_dups(example) result(results)
          implicit none
          real(kind=8) :: example(:,:)    ! The input
          real(kind=8),allocatable :: res(:,:) ! intermediate array
          real(kind=8),allocatable :: results(:,:) ! The output
          integer :: i,j,m,k   ! index
          k = 1
          m = size(example,2)
          allocate(res(2,m))
          res(1:2,1) = example(1:2,1)
          outer: do i=2,m
                  do j=1,k
                     if (dist(res(1:2,j),example(1:2,i))<1d-14) then
                    ! Found a match so start looking again
                       cycle outer
                     end if
                  end do
                  ! No match found so add it to the output
                  k = k+1
                  res(1:2,k) = example(1:2,i)
                 end do outer
          allocate(results(2,k))
          results = res(1:2,1:k)
          !write(*,*) "The result is:  ",res
         end function remove_dups

      SUBROUTINE KB07AD(COUNT,N,INDEX)
!  TO BE SORTED.

!     .. Scalar Arguments ..
      INTEGER N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M
      INTEGER MLOOP
!     ..
!     .. Local Arrays ..
      INTEGER MARK(100)

!     ..
!     .. Executable Statements ..
!  SET INDEX ARRAY TO ORIGINAL ORDER .
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
!  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED .
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)

   20 FORMAT (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ',&
      'RETURN TO CALLING PROGRAM')

      GO TO 200
!  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
!  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
   30 M = 12
!  SET UP INITIAL VALUES.
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
!  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
      IFKA = IF - IS
       IF ((IFKA+1).GT.M) GO TO 70
!********* FINAL SORTING ***
!  ( A SIMPLE BUBBLE SORT )
       IS1 = IS + 1
       DO 60 J = IS1,IF
         I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
         IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
         IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
!             *******  QUICKSORT  ********
!  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
!  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
!  HIGHEST ADDRESS.
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
!  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
!  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
!  OF X .
        K = 1
        IFK = IF
!  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
!  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS
!  NECESSARY, UNTIL THEY MEET .
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     CONTINUE
          GO TO 120

  110   CONTINUE
!  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
!  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
!  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL
!  TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED
!  INDEPENDENTLY .
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
!  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
!  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
!  IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE .
        LA = LA + 2
        GO TO 190

  170   IF (LA.LE.0) GO TO 200
!  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN

      END SUBROUTINE KB07AD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     subroutine cut_cells_find(N,x,y,x_lower,y_lower,dx,dy,&
!       mx,my,maxm,mbc,ii,jj,num_cut_cells)
!       ! this subroutine finds the index-pairs of the cut cells given a straight line of embedded barrier
!       !  input: x_0,y_0,x_e,y_e :: real(8) array of N intersecting points of barrier with grid
!       !         i_0,j_0,i_e,j_e :: integer of start and end cells' indices of barrier
!       !         x_lower,y_lower the lower left corner of grid
!       !         dx,dy are mesh sizes of grid
!       !         maxm is max(mx,my) the max number of grid cells in a row or col
!       !         mbc is number of ghost cells
!       ! output: ii the x-index of cut cells (1D array of integers)
!       !         jj the y-index of cut cells (1D array of integers)
!       !         num_cut_cells is the number of cut cells
!
!       implicit none
!       integer :: N
!       real(kind=8) :: x(N),y(N),x_lower,y_lower,dx,dy,slope_bar,theta
!       integer :: i_0,j_0,i_e,j_e,maxm, ii(N-1), jj(N-1),i,i_count,j,k
!       integer :: length,start,end, mbc,mx,my, num_cut_cells
!       integer, allocatable :: index_ord(:)
!       real(kind=8), allocatable :: intersect_cut(:,:), midpoints(:,:)
!       real(kind=8) :: xe(1-mbc:mx+mbc+1), ye(1-mbc:my+mbc+1)
!       real(kind=8) :: intersections(2*maxm,2), pi, dist_x,dist_y
!       logical :: vertical
!
!       ! barrier angle:
!       pi = 3.14159265359d0
!       vertical = .false.
!
!       ! physical grid setup
!       ! do loop to make the grid nodes 'xe' and 'ye'
!       xe(1-mbc)=x_lower - 2*dx
!       xe(1-mbc+1)=x_lower - dx
!       do i = 1, mx+mbc+1
!         xe(i) = x_lower + (i-1)*dx
!       end do
!       ye(1-mbc)=y_lower - 2*dy
!       ye(1-mbc+1)=y_lower - dy
!       do i = 1, my+mbc+1
!         ye(i) = y_lower + (i-1)*dy
!       end do
!
!       ! for each pair of neighboring intersections, find which cell it cuts
!          !find the midpoint of the intersections and see which cell they fit in
!          ! start with i_0,j_0 and walk up or down
!         allocate(midpoints(N-1,2))
!         do i=1, N-1
!           midpoints(i,1:2) = (/0.5d0*(x(i)+x(i+1)),&
!           0.5d0*(y(i)+y(i+1))/)
!         end do
!         ! write(*,*) "midpoint ", midpoints
!
!       ! put them into ii and jj
!       ! the first cell where the barrier starts is always cut:
!       do i=1-mbc,mx+mbc
!         if (midpoints(1,1) .gt. xe(i) .and. midpoints(1,1) .lt. xe(i+1)) then
!           i_0 = i
!           exit
!         end if
!       end do
!       do i = 1-mbc,my+mbc
!         if (midpoints(1,2) .gt. ye(i) .and. midpoints(1,2) .lt. ye(i+1)) then
!           j_0 = i
!           exit
!         end if
!       end do
!
!       ii(1) = i_0
!       jj(1) = j_0
!       j = 1 ! counter for cut cells
!       do i = 1, N-2
!         if (midpoints(i+1,1) .lt. xe(ii(i)) .and. midpoints(i+1,1) .gt. xe(ii(i)-1)) then
!           ii(i+1) = ii(i) - 1
!         else if (midpoints(i+1,1) .gt. xe(ii(i)+1) .and. midpoints(i+1,1) .lt. xe(ii(i)+2)) then
!           ii(i+1) = ii(i) + 1
!         else if (midpoints(i+1,1) .gt. xe(ii(i)) .and. midpoints(i+1,1) .lt. xe(ii(i)+1)) then
!           ii(i+1) = ii(i)
!         end if
!         if (midpoints(i+1,2) .lt. ye(jj(i)) .and. midpoints(i+1,2) .gt. ye(jj(i)-1)) then
!           jj(i+1) = jj(i) - 1
!         else if (midpoints(i+1,2) .gt. ye(jj(i)+1) .and. midpoints(i+1,2) .lt. ye(jj(i)+2)) then
!           jj(i+1) = jj(i) + 1
!         else if (midpoints(i+1,2) .gt. ye(jj(i)) .and. midpoints(i+1,2) .lt. ye(jj(i)+1)) then
!           jj(i+1) = jj(i)
!         end if
!       end do
!
!       ! number of cut cells (to index from 1 to it for ii and jj)
!       num_cut_cells = N-1
!
!     end subroutine cut_cells_find
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PSEUDO CODE:
! FOR EACH CUT CELL:
! get intersections between barrier and grid and the four corners of the cut cell
! find cycles
! for each edge of the cycle, form 2 h-boxes both in transverse and normal direction (extend off endpoints by length of h)
! find intersections between hbox and grid  (as previously)
! find cycles, find their indices, get h-box averages (for those that cross barrier, get appropriate avg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


   function find_centroid(cycl) result(cent)
     implicit none
     real(8):: cycl(:,:)
     real(8) :: c_x,c_y,cent(2), det, tempdet
     integer :: n,j,i
     j = 1
     n = size(cycl,2)
     c_x = 0.d0
     c_y = 0.d0
     det = 0.d0
     do i=1,n
       if (i .eq. n) then
         j = 1
       else
         j = i + 1
       end if
       tempdet = cycl(1,i) * cycl(2,j) - cycl(1,j)*cycl(2,i)
       det = det + tempdet
       c_x = c_x + (cycl(1,i)+cycl(1,j))*tempdet
       c_y = c_y + (cycl(2,i)+cycl(2,j))*tempdet
     end do
     c_x = c_x/(3*det)
     c_y = c_y/(3*det)
     cent = (/c_x,c_y/)
   end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine find_ind(cycl,rough_guess,xe,ye,mx,my,mbc,right_ind,centroid)
    ! find the index of the cell that cycle-polygon covers
    !
    ! INPUT :: cycl (list), list of N vertices (list) of the (h-box fragment) polygon, first element same with last
    !          rough_guess (list), list of integers representing index close to actual index
    !          xe, ye (array), physical grid nodes
    ! OUTPUT :: right_ind (list), list of integers representing actual index of the cycle-polygon
    implicit none
    integer :: N,rough_guess(2),mbc,mx,my,right_ind(2),rot,i,right_ind_x,right_ind_y
    real(kind=8) :: cycl(:,:),xe(-mbc:mx+mbc),ye(-mbc:my+mbc),centroid(2)
    real(kind=8) :: A, c_x,c_y, Ver3(2),Ver2(2),Ver1(2),vec2(2),vec1(2)
    real(kind=8) :: V1(2), V2(2),x_guess,y_guess

    N = size(cycl,2)
    ! # centroid coordinate calculation:
    A = area_polygon(cycl(1,:),cycl(2,:))  ! switch the dmensions atrribute of lists... for all affected parts
    c_x = 0.d0
    c_y = 0.d0
    ! determine counter or clock:
    centroid = find_centroid(cycl)
    c_x = centroid(1)
    c_y = centroid(2)
    ! rough guess and adjusting
    x_guess = xe(rough_guess(1))
    y_guess = ye(rough_guess(2))

    ! search:
    if (c_x .le. x_guess) then
        do while (c_x .le. x_guess)
            rough_guess(1) = rough_guess(1)-1
            ! print*, rough_guess(1)
            x_guess = xe(max(rough_guess(1),-mbc))
        end do
        right_ind_x = rough_guess(1)
    else if (c_x .ge. x_guess) then
        do while (c_x .ge. x_guess)
            rough_guess(1) = rough_guess(1) + 1
            x_guess = xe(min(rough_guess(1),mx+mbc))
        end do
        right_ind_x = rough_guess(1)-1
    end if
    if (c_y .le. y_guess) then
        do while (c_y .le. y_guess)
            rough_guess(2) = rough_guess(2)-1
            y_guess = ye(max(rough_guess(2),-mbc))
        end do
        right_ind_y = rough_guess(2)
    else if (c_y .ge. y_guess) then
        do while (c_y .ge. y_guess)
            rough_guess(2) = rough_guess(2) + 1
            y_guess = ye(min(rough_guess(2),my+mbc))
        end do
        right_ind_y = rough_guess(2)-1
    end if
    right_ind(1) = right_ind_x
    right_ind(2) = right_ind_y
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! function vert_ind(phys_ver,rough_guess,xe,ye,mbc,mx,my) result(right_ind)
    !     ! find index of cell in which a point lies (similar to find_ind, except this is for points and that is for polygons)
    !     !
    !     ! INPUT :: phys_ver (list), list representing a vertex coordinate
    !     !          rough_guess (list), list of integers representing index close to actual index
    !     !          xe, ye (array), physical grid nodes
    !     ! OUTPUT :: right_ind (list), list of integers representing actual index of the vertex's cell
    !     implicit none
    !     integer :: mx,my,mbc
    !     real(kind=8) :: phys_ver(2),xe(1-mbc:mx+mbc+1),ye(1-mbc:my+mbc+1)
    !     integer :: rough_guess(2), right_ind(2),right_ind_x,right_ind_y
    !     real(kind=8) :: x_guess, y_guess
    !
    !     ! rough guess and adjusting
    !     x_guess = xe(rough_guess(1))
    !     y_guess = ye(rough_guess(2))
    !
    !     ! search:
    !     if (phys_ver(1) .le. x_guess) then
    !         do while (phys_ver(1) .le. x_guess)
    !             rough_guess(1) = rough_guess(1)-1
    !             x_guess = xe(max(rough_guess(1),1))
    !         end do
    !         right_ind_x = rough_guess(1)
    !     else if (phys_ver(1) .ge. x_guess) then
    !         do while (phys_ver(1) .ge. x_guess)
    !             rough_guess(1) = rough_guess(1) + 1
    !             x_guess = xe(min(rough_guess(1),size(xe)-1))
    !         end do
    !         right_ind_x = rough_guess(1)-1
    !     end if
    !     if (phys_ver(2) .le. y_guess) then
    !         do while (phys_ver(2) .le. y_guess)
    !             rough_guess(2) = rough_guess(2)-1
    !             y_guess = ye(max(rough_guess(2),1))
    !         end do
    !         right_ind_y = rough_guess(2)
    !     else if (phys_ver(2) .ge. y_guess) then
    !         do while (phys_ver(2) .ge. y_guess)
    !             rough_guess(2) = rough_guess(2) + 1
    !             y_guess = ye(min(rough_guess(2),size(xe)-1))
    !         end do
    !         right_ind_y = rough_guess(2)-1
    !     end if
    !     right_ind(1) = right_ind_x
    !     right_ind(2) = right_ind_y
    !   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ray_intersects_seg(P,A1,B1,intersects,on_line)
    ! checks if right going ray starting from P intersects line segment A1-B1,
    ! if point P is on the line, then returns true for "on_line" bool (needed to show such point will not be
    ! counted as inside polygon)
    !
    ! INPUT :: P (list), coordinate of starting point of ray
    !         A1, B1 (list,list), coordinates of line segment's endpoints
    ! OUTPUT :: intersects (bool), true if intersects or false if not
    !           on_line (bool), true if point P is on A1-B1 or false if not
    implicit none
    integer, allocatable :: index_ord(:)
    real(kind=8) :: P(2),A1(2),B1(2), points(2,2),A(2),B(2),tol
    real(kind=8) :: m_red,m_blue
    logical :: on_line, intersects
    tol = 1d-14
    on_line = .false.
    intersects = .false.
    points(1:2,1) = A1
    points(1:2,2) = B1
    allocate(index_ord(2))
    call KB07AD(points(2,1:2),2,index_ord)
    points(1,1:2) = points(1,index_ord)
    A = points(1:2,1)
    B = points(1:2,2)
    if (abs(P(2)-A(2)).le.tol .or. abs(P(2)-B(2)).le.tol) then
        P(2) = P(2) + tol
    end if
    if (P(2) .le. A(2) .or. P(2) .gt. B(2)) then
      GOTO 99
    else if (P(1) .ge. max(A(1),B(1))) then
      GOTO 99
    else
        if (P(1) .lt. min(A(1),B(1))) then
            intersects = .true.
            GOTO 99
        else
            if (abs(A(1)-B(1)) .gt. tol) then
                m_red = (B(2)-A(2))/(B(1)-A(1))
            else
                m_red = huge(0)
            end if
            if (abs(A(1)-P(1)) .gt. tol) then
                m_blue =(P(2)-A(2))/(P(1)-A(1))
            else
                m_blue = huge(0)
            end if
            if (abs(m_blue-m_red) .lt. tol) then
                on_line = .true.
               GOTO 99
            end if
            if (m_blue .gt. m_red) then
                intersects = .true.
                GOTO 99
            else
                intersects = .false.
                GOTO 99
            end if
          end if
      end if
  99 continue
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function in_or_out(P,sides,N) result(code)
    ! finds whether point P is inside or outside of polygon made by sides
    !
    ! INPUT :: P (list), coordinate of point to check whether in polygon
    !         sides (list), collection of *N* pairs of coordinates that form a polygon
    ! OUTPUT :: code (integer), 1 if inside or -1 if outside
    implicit none
    integer :: count,code,i,N
    real(kind=8) :: P(2),A1(2),B1(2)
    real(kind=8) :: sides(N,2,2)
    logical :: intersects, on_line
    count = 0
    code = 1
    do i=1,N
      call ray_intersects_seg(P,sides(i,1:2,1),sides(i,1:2,2),intersects,on_line)
        if (on_line) then
          return
        end if
        if (intersects) then
            count = count + 1
        end if
    end do
    if (modulo(count, 2) .ne.  0) then
        GOTO 99
    else
        code = code * -1
        GOTO 99
    end if
  99 continue
    end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine is_cycle(edgs,yes_or_no)
    ! Check if set of edges form a cycle
    !
    ! INPUT :: edgs (list), collection of N pair of coordinates [[x,y],[a,b]]
    ! OUTPUT :: yes_or_no (bool), True if edgs form cycle, False if not
      implicit none
      real(real32) :: nan
      integer :: N,j,i,n1,k,j2,N2,N_uniq,N_var,k2,i3
      integer:: d1,d2
      integer ,allocatable :: indice(:)
      real(kind=8) :: edgs(:,:,:),edg(2,2),v1(2),v2(2)  ! for edgs, first index is counter, second index is (x,y), third is label for the two coordinates in an edge
      real(kind=8) :: first_edge(2,2)
      real(kind=8), allocatable :: ord_vert(:,:),unique_vert(:,:),setOfElems(:,:)
      real(kind=8), allocatable :: cycl(:,:),vertices(:,:), ord_edges(:,:,:), edgs_copy(:,:,:)
      real(kind=8) :: next_v(2),elem(2)
      logical :: yes_or_no
      ! length of list:

      N = size(edgs,3)
      d1 = size(edgs,1)
      d2 = size(edgs,2)
      allocate(edgs_copy(d1,d2,N))
      edgs_copy = edgs
      ! NaN value to indicate end of array
      nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
      ! one or two side cases automatic cancel
      if (N==1 .or. N==2) then
        yes_or_no = .false.
        return
      end if
      ! print*, "Edges:copy", edgs_copy

      !  isolate all unique vertices included in the set of edges
      allocate(vertices(2,2*N))
      vertices(1:2,(/ (j, j = 1, 2*N,2) /)) = edgs(1:2,1,(/ (i, i = 1, N) /))
      vertices(1:2,(/ (j, j = 2, 2*N,2) /)) = edgs(1:2,2,(/ (i, i = 1, N) /))
      ! intilize a null list
      unique_vert = remove_dups(vertices)

      N_uniq = size(unique_vert,2)
      N_var = N
      allocate(ord_vert(2,N_uniq+1))
      ord_vert(1:2,1) = unique_vert(1:2,1)
      i = 2
      do j=1,size(edgs_copy,3)
        edg = edgs_copy(:,:,j)
        if (dist(edg(:,1),unique_vert(:,1)).lt.1d-14.or.dist(edg(:,2),unique_vert(:,1)).lt.1d-14) then
          first_edge = edg
          n1 = j
          call IAddToList(indice,n1)
          if (dist(edg(:,1),unique_vert(:,1)).gt.1d-14) then
            ord_vert(:,i) = edg(:,1)
            i = i +1
          else if (dist(edg(:,2),unique_vert(:,1)).gt.1d-14) then
            ord_vert(:,i) = edg(:,2)
            i = i + 1
          end if
          call DelFromList_edges(edgs_copy,n1)
          exit
        end if
      enddo
      j2 = 1
      ! print*, "Edges:copy", edgs_copy
      N_var = size(edgs_copy,3)

      do while (size(edgs_copy,3) .gt. 0)
        next_v = ord_vert(:,i-1)
      do j=1,size(edgs_copy,3)
        edg = edgs_copy(1:2,1:2,j)
        ! write(*,*) "J:,",j, "EDGE CHK:", edg
        ! write(*,*) " ord vert last", ord_vert(i-1,:)
        if ((dist(edg(1:2,1),next_v).le.1d-14) .or. (dist(edg(1:2,2),next_v).le.1d-14 )) then
          n1 = j
          call IAddToList(indice,n1)
          if (dist(edg(1:2,1),next_v).gt.1d-14) then
             ord_vert(1:2,i) = edg(1:2,1)
             i = i + 1
          else if (dist(edg(1:2,2),next_v).gt.1d-14) then
            ord_vert(1:2,i) = edg(1:2,2)
            i = i + 1
          end if
          call DelFromList_edges(edgs_copy,n1)
          exit
          ! N_var = N_var - 1
        end if
        ! print*, "Edges:copy", edgs_copy

      end do
      j2 = j2 + 1
      if (j2 .gt. N_var) then
        exit
      end if
    end do
    ! write(*,*) "ORDERDE VERT x", ord_vert(1,:)
    ! write(*,*) "ordered vert y" ,ord_vert(2,:)
    ! allocate(setOfElems(2*N,2))
    ! allocate(cycl(N_uniq+1,2))
      ! if repeat of vertex in the ord vert, then there is a cycle, the first appearance to second
      yes_or_no = .false.
      k = 1
      do j = 1, i-1
        elem = ord_vert(1:2,j)
        if (j == 1) then
          call AddToList_verts(setOfElems,elem)
          ! setOfElems(j,:) = elem
          k = k + 1
          cycle
        end if
        if (vert_in_list(elem,setOfElems)) then
          call AddToList_verts(setOfElems,elem)
          yes_or_no = .true.
          exit
        else
          call AddToList_verts(setOfElems,elem)
        end if
      end do
      if (yes_or_no) then
        do j=1,size(setOfElems,2)
          if (dist(setOfElems(:,j),setOfElems(:,size(setOfElems,2))).gt.1d-14) then
            call DelFromList_nvert(setOfElems,j)
          else
            allocate(cycl(2,size(setOfElems,2)))
            cycl = setOfElems
            exit
          end if
        end do
      else
        allocate(cycl(2,size(setOfElems,2)))
        cycl = nan
      end if
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE envelope(x, y, n, vertex, nvert, iwk)
    !
    ! !  Find the vertices (in clockwise order) of a polygon enclosing
    ! !  the points (x(i), y(i), i=1, ..., n.
    !
    ! !  On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
    ! !  iwk() is an integer work array which must have dimension at least n
    ! !  in the calling program.
    !
    ! !  There is a limit of 100 vertices imposed by the dimension of array next.
    !
    ! !  Programmer: Alan Miller
    ! !  Latest revision - 12 September 1987
    ! !  Fortran 90 version - 8 August 1996
    !
    ! IMPLICIT NONE
    ! INTEGER :: n, vertex(n), nvert, iwk(n)
    ! REAL(8)    :: x(n), y(n)
    !
    ! !       Local variables
    !
    ! INTEGER :: next(100), i, i1, i2, j, jp1, jp2, i2save, i3, i2next
    ! REAL(8)    :: xmax, xmin, ymax, ymin, dista, dmax, dmin, x1, y1, dx, dy, x2, y2
    ! REAL(8)    ::           dx1, dx2, dmax1, dmax2, dy1, dy2, temp, zero = 0.0
    !
    ! IF (n < 2) RETURN
    !
    ! !  Choose the points with smallest & largest x- values as the
    ! !  first two vertices of the polygon.
    !
    ! IF (x(1) > x(n)) THEN
    !   vertex(1) = n
    !   vertex(2) = 1
    !   xmin = x(n)
    !   xmax = x(1)
    ! ELSE
    !   vertex(1) = 1
    !   vertex(2) = n
    !   xmin = x(1)
    !   xmax = x(n)
    ! END IF
    !
    ! DO i = 2, n-1
    !   temp = x(i)
    !   IF (temp < xmin) THEN
    !     vertex(1) = i
    !     xmin = temp
    !   ELSE IF (temp > xmax) THEN
    !     vertex(2) = i
    !     xmax = temp
    !   END IF
    ! END DO
    !
    ! !       Special case, xmax = xmin.
    !
    ! IF (xmax == xmin) THEN
    !   IF (y(1) > y(n)) THEN
    !     vertex(1) = n
    !     vertex(2) = 1
    !     ymin = y(n)
    !     ymax = y(1)
    !   ELSE
    !     vertex(1) = 1
    !     vertex(2) = n
    !     ymin = y(1)
    !     ymax = y(n)
    !   END IF
    !
    !   DO i = 2, n-1
    !     temp = y(i)
    !     IF (temp < ymin) THEN
    !       vertex(1) = i
    !       ymin = temp
    !     ELSE IF (temp > ymax) THEN
    !       vertex(2) = i
    !       ymax = temp
    !     END IF
    !   END DO
    !
    !   nvert = 2
    !   IF (ymax == ymin) nvert = 1
    !   RETURN
    ! END IF
    !
    ! !  Set up two initial lists of points; those points above & those below the
    ! !  line joining the first two vertices.    next(i) will hold the pointer to the
    ! !  point furthest from the line joining vertex(i) to vertex(i+1) on the left
    ! !  hand side.
    !
    ! i1 = vertex(1)
    ! i2 = vertex(2)
    ! iwk(i1) = -1
    ! iwk(i2) = -1
    ! dx = xmax - xmin
    ! y1 = y(i1)
    ! dy = y(i2) - y1
    ! dmax = zero
    ! dmin = zero
    ! next(1) = -1
    ! next(2) = -1
    !
    ! DO i = 1, n
    !   IF (i == vertex(1) .OR. i == vertex(2)) CYCLE
    !   dista = (y(i) - y1)*dx - (x(i) - xmin)*dy
    !   IF (dista > zero) THEN
    !     iwk(i1) = i
    !     i1 = i
    !     IF (dista > dmax) THEN
    !       next(1) = i
    !       dmax = dista
    !     END IF
    !   ELSE IF (dista < zero) THEN
    !     iwk(i2) = i
    !     i2 = i
    !     IF (dista < dmin) THEN
    !       next(2) = i
    !       dmin = dista
    !     END IF
    !   END IF
    ! END DO
    !
    ! !  Ends of lists are indicated by pointers to -ve positions.
    !
    ! iwk(i1) = -1
    ! iwk(i2) = -1
    ! nvert = 2
    !
    ! j = 1
    !
    ! !  Start of main process.
    !
    ! !  Introduce new vertex between vertices j & j+1, if one has been found.
    ! !  Otherwise increase j.   Exit if no more vertices.
    !
    ! 40 IF (next(j) < 0) THEN
    ! IF (j == nvert) RETURN
    ! j = j + 1
    ! GO TO 40
    ! END IF
    !
    ! jp1 = j + 1
    ! DO i = nvert, jp1, -1
    !   vertex(i+1) = vertex(i)
    !   next(i+1) = next(i)
    ! END DO
    ! jp2 = jp1 + 1
    ! nvert = nvert + 1
    ! IF (jp2 > nvert) jp2 = 1
    ! i1 = vertex(j)
    ! i2 = next(j)
    ! i3 = vertex(jp2)
    ! vertex(jp1) = i2
    !
    ! !  Process the list of points associated with vertex j.   New list at vertex j
    ! !  consists of those points to the left of the line joining it to the new
    ! !  vertex (j+1).   Similarly for the list at the new vertex.
    ! !  Points on or to the right of these lines are dropped.
    !
    ! x1 = x(i1)
    ! x2 = x(i2)
    ! y1 = y(i1)
    ! y2 = y(i2)
    ! dx1 = x2 - x1
    ! dx2 = x(i3) - x2
    ! dy1 = y2 - y1
    ! dy2 = y(i3) - y2
    ! DMAX1 = zero
    ! dmax2 = zero
    ! next(j) = -1
    ! next(jp1) = -1
    ! i2save = i2
    ! i2next = iwk(i2)
    ! i = iwk(i1)
    ! iwk(i1) = -1
    ! iwk(i2) = -1
    !
    ! 60 IF (i /= i2save) THEN
    !   dista = (y(i) - y1)*dx1 - (x(i) - x1)*dy1
    !   IF (dista > zero) THEN
    !     iwk(i1) = i
    !     i1 = i
    !     IF (dista > DMAX1) THEN
    !       next(j) = i
    !       DMAX1 = dista
    !     END IF
    !   ELSE
    !     dista = (y(i) - y2)*dx2 - (x(i) - x2)*dy2
    !     IF (dista > zero) THEN
    !       iwk(i2) = i
    !       i2 = i
    !       IF (dista > dmax2) THEN
    !         next(jp1) = i
    !         dmax2 = dista
    !       END IF
    !     END IF
    !   END IF
    !   i = iwk(i)
    ! ELSE
    !   i = i2next
    ! END IF
    !
    ! !  Get next point from old list at vertex j.
    !
    ! IF (i > 0) GO TO 60
    !
    ! !  End lists with -ve values.
    !
    ! iwk(i1) = -1
    ! iwk(i2) = -1
    !
    ! GO TO 40
    ! END SUBROUTINE envelope

    subroutine line_crossing(m1,b1,m2,b2,x_star,y_star,code)
      implicit none
      real(8) :: m1,m2,b1,b2,x_star,y_star
      integer :: code

      ! input: slopes m1, m2, and y-ints b1, b2
      ! output: intersection x_star,y_star, and integer code: 1, found intersection, 2, parallel

      if (abs(m1-m2).le.1d-14) then
        code = 2
        x_star =0.d0
        y_star =0.d0
        return
      else
        x_star = (b2-b1)/(m1-m2)
        y_star = m1*x_star + b1
        code =1
      end if
    end subroutine

     subroutine edge_in_region(edge,region,N,true)
       implicit none
       integer :: N,i
       real(8) :: edge(2,2), region(N,2,2)
       logical :: true

       do i = 1,N
         if (dist(edge(1:2,1),region(i,1:2,1)).le.1d-14 .and. dist(edge(1:2,2),region(i,1:2,2)).le.1d-14) then
           true = .true.
           exit
         else
           true = .false.

         end if
       end do
     end subroutine

     function vert_in_edge(vert,edge) result(true)
       implicit none
       real(8) :: vert(2), edge(2,2)
       logical :: true

       if (dist(vert,edge(1:2,1)).le.1d-14 .or. dist(vert,edge(1:2,2)).le.1d-14) then
         true = .true.
       else
         true = .false.
       end if
     end function

     function ind_in_inds(ind,inds) result(truea)
       implicit none
       integer :: ind(2),inds(:,:), i,N
       logical :: truea
       N = size(inds,2)
       truea=.false.
       do i=1,N
         if (ind(1)==inds(1,i) .and. ind(2)==inds(2,i)) then
           truea = .true.
         end if
       end do
     end function


    function orientation(p,q,r) result(o)
      ! input : verts p,q,r
      ! output : integer o (0 if colinear, 1 if CW, 2 if CCW)
      implicit none
      real(8) :: p(2),q(2),r(2)
      integer :: o
      ! local ::
      real(8) :: val

      val = (q(2)-p(2))*(r(1)-q(1)) - (q(1)-p(1))*(r(2)-q(2))
      if (val .gt. 0.d0) then
        o = 1
      else if (val.lt.0.d0) then
        o = 2
      else
        o = 0
      end if
    end function

    function onSegment(p,q,r) result(yes)
      ! input : verts p,q,r
      ! output:  logical yes (true if q lies on pr, false if not)
      implicit none
      real(8) :: p(2),q(2),r(2)
      logical :: yes
      if (q(1).le.max(p(1),r(1)) .and. q(1).ge.min(p(1),r(1)) .and. q(2).le.max(p(2),r(2)) .and. q(2).ge.min(p(2),r(2))) then
        yes = .true.
      else
        yes = .false.
      end if
    end function

    function doIntersect(p1,q1,p2,q2) result(yes)
      ! input : verts p1,q1,p2,q2, the same numbered ones form a segment
      ! output : logical yes, true if two line segs intersect and false if not
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2)
      logical :: yes
      ! local ::
      integer :: o1,o2,o3,o4

      o1 = orientation(p1,q1,p2)
      o2 = orientation(p1,q1,q2)
      o3 = orientation(p2,q2,p1)
      o4 = orientation(p2,q2,q1)

      yes = .false. ! default
      if ((o1 /= o2) .and. (o3 /= o4)) then
        yes = .true.
      end if


      if ((o1 == 0) .and. onSegment(p1,p2,q1)) then
        yes = .true.
      else if ((o2 == 0) .and. onSegment(p1, q2, q1)) then
       yes = .true.
     else if ((o3 == 0) .and. onSegment(p2, p1, q2)) then
       yes = .true.
     else if ((o4 == 0) .and. onSegment(p2, q1, q2)) then
       yes = .true.
     end if

    end function

    subroutine lineseg_intersections(p1,q1,p2,q2,yes,p_star)
      ! input: verts p1,q1,p2,q2
      ! output: logical yes if intersect , and vert p_star, the intersection
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2),p_star(2)
      logical :: yes
      ! local
      real(8) :: m1,b1,m2,b2,x_star,y_star
      integer :: code

      yes = doIntersect(p1,q1,p2,q2)
      if (abs(min(p1(1),q1(1))-min(p2(1),q2(1))).le.1d-14 .or. abs(max(p1(1),q1(1)) - max(p2(1),q2(1))).le.1d-14) then
        yes = .false.
      else if (abs(min(p1(2),q1(2))-min(p2(2),q2(2))).le.1d-14 .or. abs(max(p1(2),q1(2))-max(p2(2),q2(2))).le.1d-14) then
        yes = .false.
      end if
      if (yes) then
        if (abs(q1(1)-p1(1)).le.1d-14) then
          x_star = p1(1)
          y_star = (q2(2)-p2(2))/(q2(1)-p2(1)) * x_star + (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
        else if (abs(q2(1)-p2(1)).le.1d-14) then
          x_star = p2(1)
          y_star = (q1(2)-p1(2))/(q1(1)-p1(1)) * x_star + (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
        else
          m1 = (q1(2)-p1(2))/(q1(1)-p1(1))
          m2 = (q2(2)-p2(2))/(q2(1)-p2(1))
          b1 = (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
          b2 = (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
          call line_crossing(m1,b1,m2,b2,x_star,y_star,code)
        end if
      else if (.not. yes) then
        return
      end if
    end subroutine

    subroutine lineseg_intersections_endpoint(p1,q1,p2,q2,yes,p_star)
      ! input: verts p1,q1,p2,q2
      ! output: logical yes if intersect , and vert p_star, the intersection
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2),p_star(2)
      logical :: yes
      ! local
      real(8) :: m1,b1,m2,b2,x_star,y_star
      integer :: code

      yes = doIntersect(p1,q1,p2,q2)
      if (yes) then
        if (abs(q1(1)-p1(1)).le.1d-14) then
          x_star = p1(1)
          y_star = (q2(2)-p2(2))/(q2(1)-p2(1)) * x_star + (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
        else if (abs(q2(1)-p2(1)).le.1d-14) then
          x_star = p2(1)
          y_star = (q1(2)-p1(2))/(q1(1)-p1(1)) * x_star + (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
        else
          m1 = (q1(2)-p1(2))/(q1(1)-p1(1))
          m2 = (q2(2)-p2(2))/(q2(1)-p2(1))
          b1 = (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
          b2 = (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
          call line_crossing(m1,b1,m2,b2,x_star,y_star,code)
        end if
      else if (.not. yes) then
        return
      end if
    end subroutine


    ! subroutine find_cycles(start_edge,edges_list,N,hbox_area,frag_ratio,frag_inds,ii,jj, &
    !                              xe,ye,mx,my,mbc,orient,centroid,frag_cycle,len_cycle)
    !   implicit none
    !   ! given list of edges, finds a cycle that has start_edge as its first edge, in one orientation (sign orient)
    !   integer :: N ,frag_inds(2),mx,my,mbc,ii,jj,orient,len_cycle
    !   real(8) :: start_edge(2,2),edges_list(N,2,2), hbox_area, frag_ratio
    !   real(8) :: xe(1-mbc:mx+mbc+1),ye(1-mbc:my+mbc+1),centroid(2),frag_cycle(10,2)
    !   ! local :
    !   integer :: i,j,j2,j3,j4,j5,right_inds(2),rot,n1,exclud(N),e1
    !   real(8) :: cycle_edges(10,2,2),cycle_verts(11,2),edge(2,2),v1(2),v2(2),x1,x2,y1,y2
    !   real(8) :: last_elt(2),v_last1(2),v_last2(2),cand(2),Ver3(2),Ver2(2),Ver1(2),vec2(2),vec1(2)
    !   real(8) :: z,area,c_x,c_y
    !   logical :: yes,true1,true2,yes2
    !
    !     j2 = 1
    !     j4 = 1
    !     j5 = 1
    !     exclud=0
    !     e1 =1
    !     cycle_edges(j2,:,:) = start_edge
    !     cycle_verts(j4,:) = start_edge(1:2,1)
    !     cycle_verts(j4+1,:) = start_edge(1:2,2)
    !     j2 = j2 + 1 ! cycle edge counter
    !     j4 = j4 + 2 ! cycle vertices counter
    !     yes = .false.
    !     do while (.not. yes)
    !       ! write(*,*) "excluded ", exclud
    !       do j3 = 1,N
    !       !   call int_in_list(j3,exclud,N,yes2)
    !       !   if (yes2) then
    !       !     cycle
    !       !   end if
    !         edge = edges_list(j3,:,:)
    !         v1 = edge(1:2,1)
    !         v2 = edge(1:2,2)
    !         x1 = v1(1)
    !         x2 = v2(1)
    !         y1 = v1(2)
    !         y2 = v2(2)
    !         last_elt = cycle_verts(j4-1,:)
    !         ! v_last1 = last_elt(1:2,1)
    !         ! v_last2 = last_elt(1:2,2)
    !         if ((dist(v1,last_elt).le.1d-8 .or. dist(v2,last_elt).le.1d-8)) then ! .or. &
    !             ! (dist(v1,v).le.1d-8 .or. dist(v2,v_last1).le.1d-8)) then
    !           if (dist(v1,last_elt).le.1d-8) then
    !             cand = v2
    !           else if (dist(v2,last_elt).le.1d-8) then
    !             cand = v1
    !           ! else if (dist(v1,v_last1).le.1d-8) then
    !           !   cand = v2
    !           ! else if (dist(v2,v_last1).le.1d-8) then
    !           !   cand = v1
    !           end if
    !           Ver3 = cand
    !           Ver2 = cycle_verts(j4-1,:)
    !           Ver1 = cycle_verts(j4-2,:)
    !           vec2 = Ver3 - Ver2
    !           vec1 = Ver2 - Ver1
    !          ! determine orientation of rotation
    !          z = vec2(1) * vec1(2) - vec1(1)*vec2(2)
    !          write(*,*) "cand:",cand,j3
    !          if (z*orient .gt. 0 .and. abs(z) .gt. 1d-12) then
    !            ! call vert_in_edge(cycle_verts(j4-1,:),edge,true1)
    !            call edge_in_region(edge,cycle_edges(1:j2-1,:,:),j2-1,true2)
    !            if (.not. true2) then
    !              cycle_verts(j4,:)=cand
    !              cycle_edges(j2,:,:)=edge
    !              write(*,*) "EDGE added: ", j3
    !              call is_cycle(cycle_edges(1:j2,:,:),yes,cycle_verts(1:j4,:),j2,n1)
    !              j4 = j4 + 1
    !              j2 = j2 + 1
    !              ! exclud(e1) = j3
    !              ! e1 = e1 + 1
    !              exit
    !            end if
    !          end if
    !          continue
    !        end if
    !        end do
    !        write(*,*) "IN PROG: x", cycle_verts(1:j4-1,1)
    !        write(*,*) "IN PROG: y", cycle_verts(1:j4-1,2)
    !      end do
    !      area = area_polygon(cycle_verts(:,1:j4-1))
    !      call find_ind(cycle_verts(1:j4-1,:),j4-1,(/ii,jj/),xe,ye,mx,my,mbc,frag_inds,centroid)
    !      frag_ratio = area/hbox_area
    !      len_cycle = j4-1
    !      frag_cycle(1:j4-1,:) = cycle_verts(1:j4-1,:)
    !    end subroutine

      function same_side_as_small_cell(up,centroid,m,b) result(code)
        implicit none
        logical :: up
        real(8) :: centroid(2),m,b
        integer :: code

        if (up) then
          if ((centroid(2) .gt. m*centroid(1) + b)) then
            code = 1
          else
            code = -1
          end if
        else
          if ((centroid(2) .lt. m*centroid(1) + b)) then
            code = 1
          else
            code = -1
          end if
        end if
      end function

      function reflect(P,m,b) result(P_prime)
        ! point P to be reflected; mx+b line of equation to be reflected across; P_prime the reflection
        implicit none
        real(8) :: P(2),m,b,P_prime(2)
        real(8) :: m_s,t,L(2),delx,dely
        ! check if point P is actually on the line, then return just P back
        if (abs(P(2)-m*P(1)-b).le.1d-14) then
          P_prime = P
          return
        end if
        ! all other cases:
        m_s = -1/m
        t= P(2) - m_s*P(1)
        L(1) = (b-t)/(m_s-m)
        L(2) = L(1)*m_s + t
        delx = (L(1)-P(1))
        dely = (L(2)-P(2))
        P_prime(1) = P(1) + 2*delx
        P_prime(2) = P(2) + 2*dely
      end function

      function colinear_check(P1,P2,P3) result(yes)
        ! check if points P1,P2,P3 are colinear, returns yes=.true. if so, and if not, yes=.false.
        implicit none
        real(8) :: P1(2),P2(2),P3(2)
        logical :: yes
        integer :: i
        real(8) :: A, x1,x2,x3,y1,y2,y3

        yes = .false.
        x1 = P1(1)
        y1 = P1(2)
        x2 = P2(1)
        y2 = P2(2)
        x3 = P3(1)
        y3 = P3(2)

        A = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
        if (abs(A).le.1d-14) then
          yes = .true.
        end if
      end function

      subroutine int_in_list(int,list,N,yes)
        implicit none
        integer :: int, N, list(N), i
        logical :: yes
        yes = .false.
        do i = 1,N
          if (int == list(i)) then
            yes = .true.
          end if
        end do
      end subroutine

      function vert_in_list(vert, list) result(yes_or_no)
        implicit none
        real(8) :: vert(2)
        real(8) :: list(:,:)
        integer :: N ,i
        logical :: yes_or_no
        N = size(list,2)
        yes_or_no = .false.
        do i = 1,N
          if (dist(vert,list(:,i)).lt.1d-7) then
            yes_or_no = .true.
            return
          end if
        end do
        return
      end function

end module aux_module_hbox

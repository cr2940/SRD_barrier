program check_aux_SRD

use aux_module_SRD
use hbox_layer

implicit none
integer, parameter :: mx=150,my=150,mbc=2
real(8) :: x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,intersections(2,mx*4)
real(8) :: xupper,yupper,x_1,y_1,x_2,y_2,intersections2(2,mx*4)
integer :: ii(mx*4),jj(mx*4),N,type_supper(mx*4),type_sunder(mx*4),m,m2
integer :: ii2(mx*4),jj2(mx*4),N2,type_supper2(mx*4),type_sunder2(mx*4)
integer :: all_undercells_i(mx*4), all_undercells_j(mx*4)
integer :: unS_cells_i(mx*4), unS_cells_j(mx*4)
integer :: upS_cells_i(mx*4), upS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
integer :: unS_cells_i2(mx*4), unS_cells_j2(mx*4)
integer :: upS_cells_i2(mx*4), upS_cells_j2(mx*4) ! THESE will tell you who the neighbor is for small cells
integer :: all_uppercells_i(mx*4), all_uppercells_j(mx*4),k_count_un,k_count_up
integer :: N_ij_un(1-mbc:mx+mbc,1-mbc:my+mbc),N_ij_up(1-mbc:mx+mbc,1-mbc:my+mbc)
real(8) :: lengths_supper(5,4*mx), lengths_sunder(5,4*mx)
real(8) :: lengths_supper2(5,4*mx), lengths_sunder2(5,4*mx)
real(8) :: area_supper(4*mx), area_sunder(4*mx),area_supper2(4*mx), area_sunder2(4*mx)
real(8) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
integer :: num_frags_d(mx*4),index_frags_d(2,4,mx*4),num_frags_d2(mx*4),index_frags_d2(2,4,mx*4)
integer :: num_frags_u(mx*4),index_frags_u(2,4,mx*4),num_frags_u2(mx*4),index_frags_u2(2,4,mx*4)
real(8) :: hbox_areas_u(mx*4),area_frags_u(4,mx*4),hbox_areas_u2(mx*4),hbox_areas_d2(mx*4)
real(8) :: hbox_areas_d(mx*4),area_frags_d(4,mx*4),area_frags_d2(4,mx*4),area_frags_u2(4,mx*4)
x_0 = 0.0d0
y_0 = .72d0!1.d0 0.719
x_e = 0.5d0 !0.995d0
y_e = 0.412d0!0.005d0
x_1 = x_e
x_2 = 1.d0
y_1 = y_e
y_2 = y_0
dx = 0.00666666d0
dy = 0.00666666d0
xlower=0.d0
xupper = 1.d0
ylower=0.d0
yupper = 1.d0
call cut_cells_find(x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,mx,my,ii,jj,intersections,N)
call cut_cells_find(x_1,y_1,x_2,y_2,dx,dy,xlower,ylower,mx,my,ii2,jj2,intersections2,N2)

call small_cells_geom(ii(1:N),jj(1:N),intersections(:,1:N+1),mx,my,dx,dy,xlower,ylower,N,&
  type_supper,type_sunder,lengths_supper,lengths_sunder,area_supper,area_sunder)
call small_cells_geom(ii2(1:N2),jj2(1:N2),intersections2(:,1:N2+1),mx,my,dx,dy,xlower,ylower,N2,&
  type_supper2,type_sunder2,lengths_supper2,lengths_sunder2,area_supper2,area_sunder2)

call SRD_undercells(N,ii(1:N),jj(1:N),mx,my,area_sunder(1:N),x_0,y_0,x_e,y_e,dx,dy,&
  unS_cells_i,unS_cells_j,N_ij_un,all_undercells_i,all_undercells_j,k_count_un)
call SRD_undercells(N2,ii2(1:N2),jj2(1:N2),mx,my,area_sunder2(1:N2),x_1,y_1,x_2,y_2,dx,dy,&
  unS_cells_i2,unS_cells_j2,N_ij_un,all_undercells_i,all_undercells_j,k_count_un)

call SRD_uppercells(N,ii(1:N),jj(1:N),mx,my,area_supper(1:N),x_0,y_0,x_e,y_e,dx,dy,&
    upS_cells_i,upS_cells_j,N_ij_up,all_uppercells_i,all_uppercells_j,k_count_up)
call SRD_uppercells(N2,ii2(1:N2),jj2(1:N2),mx,my,area_supper2(1:N2),x_1,y_1,x_2,y_2,dx,dy,&
    upS_cells_i2,upS_cells_j2,N_ij_up,all_uppercells_i,all_uppercells_j,k_count_up)

call area_cells(ii(1:N),jj(1:N),mx,my,mbc,area_sunder(1:N),area_supper(1:N),up_area_ij,un_area_ij,.false.)
call area_cells(ii2(1:N2),jj2(1:N2),mx,my,mbc,area_sunder2(1:N2),area_supper2(1:N2),up_area_ij,un_area_ij,.true.)
! try writing the data into a text file:
open (unit=1,file = "./small_cells_data.txt")
write (1,*) N
write (1,*) N2
write (1,*) ii(1:N)
write (1,*) jj(1:N)
write (1,*) ii2(1:N2)
write (1,*) jj2(1:N2)
write (1,*) intersections(:,1:N+1)
write (1,*) intersections2(:,1:N2+1)
write (1,*) type_supper(1:N)
write (1,*) type_sunder(1:N)
write (1,*) type_supper2(1:N2)
write (1,*) type_sunder2(1:N2)
write (1,*) lengths_supper(:,1:N)
write (1,*) lengths_sunder(:,1:N)
write (1,*) lengths_supper2(:,1:N2)
write (1,*) lengths_sunder2(:,1:N2)
write (1,*) unS_cells_i(1:N)
write (1,*) unS_cells_j(1:N)
write (1,*) unS_cells_i2(1:N2)
write (1,*) unS_cells_j2(1:N2)
write (1,*) upS_cells_i(1:N)
write (1,*) upS_cells_j(1:N)
write (1,*) upS_cells_i2(1:N2)
write (1,*) upS_cells_j2(1:N2)
write (1,*) area_supper(1:N)
write (1,*) area_sunder(1:N)
write (1,*) area_supper2(1:N2)
write (1,*) area_sunder2(1:N2)
write (1,*) xlower
write (1,*) ylower
write (1,*) xupper
write (1,*) yupper
write (1,*) x_0
write (1,*) y_0
write (1,*) x_1
write (1,*) y_1
write (1,*) x_e
write (1,*) y_e
write (1,*) x_2
write (1,*) y_2
write (1,*) up_area_ij
write (1,*) un_area_ij
ENDFILE 1
REWIND 1
close (1,status="keep")

call down_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii(1),jj(1),x_0,y_0,ii(N),&
   jj(N),x_e,y_e,m,hbox_areas_d,num_frags_d,area_frags_d,index_frags_d,1)
call down_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii2(1),jj2(1),x_1,y_1,ii2(N2),&
  jj2(N2),x_2,y_2,m2,hbox_areas_d2,num_frags_d2,area_frags_d2,index_frags_d2,1)
call up_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii(1),jj(1),x_0,y_0,ii(N),&
   jj(N),x_e,y_e,m,hbox_areas_u,num_frags_u,area_frags_u,index_frags_u,1)
call up_hboxes(xlower,ylower,xupper,yupper,mbc,mx,my,ii2(1),jj2(1),x_1,y_1,ii2(N2),&
  jj2(N2),x_2,y_2,m2,hbox_areas_u2,num_frags_u2,area_frags_u2,index_frags_u2,1)
open (unit=2,file="./hbox_data.txt")
write(2,*) m
write(2,*) m2
write(2,*) num_frags_d
write(2,*) num_frags_d2
write(2,*) num_frags_u
write(2,*) num_frags_u2
write(2,*) index_frags_d
write(2,*) index_frags_d2
write(2,*) index_frags_u
write(2,*) index_frags_u2
write(2,*) hbox_areas_d
write(2,*) hbox_areas_d2
write(2,*) hbox_areas_u
write(2,*) hbox_areas_u2
write(2,*) area_frags_u
write(2,*) area_frags_u2
write(2,*) area_frags_d
write(2,*) area_frags_d2
ENDFILE 2
REWIND 2
close(2,status="keep")
! read(1,*) M
! read(1,*) ii(1:N)
print*, "The number of cut cells is: ", N
print*, "Cut cells i-index is given by:", ii(1:N)
print*, "Cut cells j-index is given by:", jj(1:N)
print*, "the intersections x are :", intersections(1,1:N+1)
print*, "the intersections y are :", intersections(2,1:N+1)
print*, "TYPES:UP",type_supper(1:N)
print*, "TYPES:UN",type_sunder(1:N)
print*, "Length ratio to dx are: UP", lengths_supper(:,1:N)
print*, "Length ratio to dx are: UN", lengths_sunder(:,1:N)
print*, "AREAs : UP", up_area_ij(ii(3),jj(3))
print*, "AREAs ; UN", un_area_ij(ii(3),jj(3))
print*, "TOTAL AREA :", area_supper(1:N) + area_sunder(1:N)
end program

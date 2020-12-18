! ========================================================================================
! 2D H-box Module - contains functions to form h-box arrays both above/right or below/left
! of barrier.
! =========================================================================================

module hbox2D

      implicit none

contains

         subroutine down_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,layer,&
                    qold,aux,meqn,num_aux,mbc,mx,my,maxm,dx,dy,q_hbox,aux_hbox,alpha)
           ! forms layer of h-boxes below the barrier (or left of vertical barrier)
           ! Input: type(point_info) begin_pt = where barrier begins, end_pt = where barrier ends
           !        integer layer = 1 for first layer of hbox or 2 for second layer
           !        qold = initial q values
           ! Output: real(3,:) q_hbox = hbox values
         implicit none
         ! declare variables
         real(kind=8) :: dx, xc(1-mbc:mx+mbc), dy,yc(1-mbc:mx+mbc)
         real(kind=8) :: x_0,y_0,x_e,y_e,pi,xlower,ylower
         real(kind=8) :: slope_bar, theta, alpha, dist_x,dist_y
         real(kind=8) :: xe(1-mbc:mx+mbc+1), ye(1-mbc:my+mbc+1),x_1,y_1
         real(kind=8) :: qold(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
         real(kind=8) :: aux(num_aux,1-mbc:mx+mbc,1-mbc:my+mbc)
         real(kind=8), allocatable :: intersections(:,:), ver_hbox(:,:)
         real(kind=8), allocatable :: ver_hbox2(:,:),we(:,:),wc(:,:)
         real(kind=8),allocatable :: we_1(:,:),wc_1(:,:),hbox_cen(:,:)
         integer :: i_0_p, i_0_t, j_0_p, j_0_t, m,j_1, j_2, j0, i_2,i_1
         integer :: j, n, i_0,j_0,i_e,j_e,i,layer,start,end,length
         integer :: mbc,mx,my,maxm,meqn,l1,l2,l3,l4,num_aux
         integer,allocatable :: rough_ind(:,:),no_edges(:),odd(:),even(:)
         integer,allocatable :: index_ord(:)
         logical :: vertical
         real(kind=8), intent(out) :: q_hbox(3,1-mbc:3*maxm), aux_hbox(num_aux,1-mbc:3*maxm)
         integer :: m1,m2

         ! pi:
         pi = 3.14159265359d0
         ! domain :
         xlower =0.d0
         ylower = 0.d0

         ! physical grid setup
         ! do loop to make the grid nodes 'xe' and 'ye'
         xe(1-mbc)=xlower - 2*dx
         xe(1-mbc+1)=xlower - dx
         do i = 1, mx+mbc+1
           xe(i) = xlower + (i-1)*dx
         end do
         ye(1-mbc)=ylower - 2*dy
         ye(1-mbc+1)=ylower - dx
         do i = 1, my+mbc+1
           ye(i) = ylower + (i-1)*dy
         end do
         xc(1-mbc:mx+mbc) = xe(1-mbc:mx+mbc) + dx/2 ! center nodes
         yc(1-mbc:my+mbc) = ye(1-mbc:my+mbc) + dy/2

         ! indexing needed later
         i_0_p = i_0
         i_0_t = i_0_p
         j_0_p = j_0
         j_0_t = j_0_p

          ! barrier slope
         vertical = .false.

         if (abs(x_e - x_0).gt.1d-7) then
           slope_bar = (y_e-y_0)/(x_e-x_0)
           theta = atan(slope_bar)
         else
           slope_bar = HUGE(1.D0)
           vertical = .true.
         end if
   !================================================
   !  ------------------------------------  horizontal case
   !================================================
         if (abs(slope_bar).lt.1d-7) then
            if (j_0 /= j_e) then
              write(*,*) "The wall is not horizontal"
            end if

            ! distance away from a physical edge
            alpha = abs(y_0 - ye(j_0+1))/dy
            if (alpha.lt.1d-7) then  ! if close enough to edge
              alpha = 0.d0
            end if
         ! down hboxes
         if (j_0 == 1) then
           write(*,*) "Wall is too high, h-box will be off the grid"
           return
         end if
         if (j_0 == 2) then
           write(*,*) "Can only make one layer of h-box"
           if (layer==2) then
             return
           end if
         end if
         j_2 = j_0 - 2
         j_1 = j_0 - 1
         if (layer == 1) then
           j0 = i_0
           ! do j = 1,m
         q_hbox(1:3,1-mbc:mx+mbc) = alpha*qold(1:3,:,j_0)+(1-alpha)*qold(1:3,:,j_1)
         aux_hbox(1,1-mbc:mx+mbc) = alpha*aux(1,:,j_0) + (1-alpha)*aux(1,:,j_1)
              ! j0 = j0 + 1
           ! end do
         end if
         if (layer == 2) then
          j0 = i_0
          ! do j = 1,m
         q_hbox(1:3,1-mbc:mx+mbc) = alpha*qold(1:3,:,j_1)+(1-alpha)*qold(1:3,:,j_2)
         aux_hbox(1,1-mbc:mx+mbc) = alpha*aux(1,:,j_1) + (1-alpha)*aux(1,:,j_2)
            ! j0 = j0 + 1
          ! end do
         end if
         return

         end if  ! end of if flat barrier case
   !=========================================================
   !                     |
   !                     |      vertical case
   !                     |
   !=========================================================
         if (vertical) then
           if (i_0 /= i_e) then
             write(*,*) "The wall is not vertical"
             return
           end if

            ! distance away from physical edge
            alpha = abs(x_0 - xe(i_0+1))/dx
            if (alpha.lt.1d-7) then
              alpha = 0.d0
            end if

            ! the left h-boxes
            if (i_0 == 1) then
              write(*,*) "Wall is too left, h-box will be off the grid"
              return
            end if
            if (i_0 == 2) then
              write(*,*) "Can only make one layer of h-box"
              if (layer==2) then
                return
              end if
            end if
            i_2 = i_0 - 2
            i_1 = i_0 - 1
            if (layer == 1) then
              j0 = j_0
              ! do j = 1,m
         q_hbox(1:3,1-mbc:mx+mbc) = alpha*qold(1:3,i_0,:)+(1-alpha)*qold(1:3,i_1,:)
         aux_hbox(1,1-mbc:mx+mbc) = alpha*aux(1,i_0,:) + (1-alpha)*aux(1,i_1,:)
                ! j0 = j0 + 1
              ! end do
            end if
            if (layer == 2) then
              j0 = j_0
              ! do j = 1,m
         q_hbox(1:3,1-mbc:mx+mbc) = alpha*qold(1:3,i_1,:)+(1-alpha)*qold(1:3,i_2,:)
         aux_hbox(1,1-mbc:mx+mbc) = alpha*aux(1,i_1,:) + (1-alpha)*aux(1,i_2,:)
                ! j0 = j0 + 1
              ! end do
            end if
            return
          end if ! end of vertical wall case

         end subroutine down_hboxes

       subroutine up_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,layer,&
              qold,aux,meqn,num_aux,mbc,mx,my,maxm,dx,dy,q_hbox,aux_hbox,alpha)
         ! forms layer of h-boxes above the barrier (or right of vertical barrier)
         ! Input: type(point_info) begin_pt = where barrier begins, end_pt = where barrier ends
         !        integer layer = 1 for first layer of hbox or 2 for second layer
         !        qold = initial q values
         ! Output: real(3,:) q_hbox = hbox values
         implicit none
         ! declare variables
         real(kind=8) :: dx, xc(1-mbc:mx+mbc), dy,yc(1-mbc:mx+mbc)
         real(kind=8) :: slope_bar, theta, alpha, pi,xlower,ylower
         real(kind=8) :: xe(1-mbc:mx+mbc+1), ye(1-mbc:my+mbc+1)
         real(kind=8) :: qold(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
         real(kind=8) :: aux(num_aux,1-mbc:mx+mbc,1-mbc:my+mbc)
         real(kind=8) :: x_0,y_0,x_e,y_e
         real(kind=8), allocatable :: intersections(:,:), ver_hbox(:,:)
         real(kind=8), allocatable :: ver_hbox2(:,:)
         integer :: i_0_p, i_0_t, j_0_p, j_0_t, m,j_1, j_2, j0, i_2,i_1
         integer :: j, n, i_0,j_0,i_e,j_e,i,layer,mbc,mx,my,maxm,meqn,i_3
         integer,allocatable :: rough_ind(:,:),no_edges(:),odd(:),even(:)
         logical :: vertical
         real(kind=8), intent(out) :: q_hbox(3,1-mbc:3*maxm), aux_hbox(num_aux,1-mbc:3*maxm)
         integer :: m1,m2,num_aux

         ! pi:
         pi = 3.14159265359d0
         ! domain:
         xlower = 0.d0
         ylower = 0.d0
         ! physical grid setup
         ! do loop to make the grid nodes 'xe' and 'ye'
         xe(1-mbc) = xlower - 2*dx
         xe(1-mbc+1) = xlower - dx
         ye(1-mbc) = ylower -2*dy
         ye(1-mbc+1) = ylower - dy
         do i = 1, mx+mbc+1
           xe(i) = xlower + (i-1)*dx
         end do
         do i = 1, my+mbc+1
           ye(i) = ylower + (i-1)*dy
         end do
         xc(1-mbc:mx+mbc) = xe(1-mbc:mx+mbc) + dx/2 ! center nodes
         yc(1-mbc:my+mbc) = ye(1-mbc:my+mbc) + dy/2

         ! indexing needed later
         i_0_p = i_0
         i_0_t = i_0_p
         j_0_p = j_0
         j_0_t = j_0_p

          ! barrier slope
         vertical = .false.

         if (abs(x_e - x_0).gt.1d-7) then
           slope_bar = (y_e-y_0)/(x_e-x_0)
           theta = atan(slope_bar)
         else
           slope_bar = HUGE(1.D0)
           vertical = .true.
         end if
   !================================================
   !  ------------------------------------ horizontal case
   !================================================
         if (abs(slope_bar).lt.1d-7) then
            if (j_0 /= j_e) then
              write(*,*) "The wall is not horizontal"
            end if

            ! distance away from a physical edge
            alpha = abs(y_0 - ye(j_0+1))/dy
            if (alpha.lt.1d-7) then  ! if close enough to edge
              alpha = 0.d0
            end if

         ! Up hboxes
         if (j_0 == size(ye)) then
           write(*,*) "Wall is too high, h-box will be off the grid"
           return
         end if
         if (j_0 == size(ye)-1) then
           write(*,*) "Can only make one layer of h-box"
           if (layer==2) then
             return
           end if
         end if
         j_2 = j_0 + 2
         j_1 = j_0 + 1
         if (layer == 1) then
           j0 = i_0
           ! do j = 1,m
       q_hbox(1:3,1-mbc:mx+mbc) = (1-alpha)*qold(1:3,:,j_1)+alpha*qold(1:3,:,j_2)
       aux_hbox(1,1-mbc:mx+mbc) = (1-alpha)*aux(1,:,j_1)+alpha*aux(1,:,j_2)
              ! j0 = j0 + 1
           ! end do
         end if
         if (layer == 2) then
          j0 = i_0
          ! do j = 1,m
       q_hbox(1:3,1-mbc:mx+mbc) = (1-alpha)*qold(1:3,:,j_2)+alpha*qold(1:3,:,j_2+1)
       aux_hbox(1,1-mbc:mx+mbc) = (1-alpha)*aux(1,:,j_2)+alpha*aux(1,:,j_2+1)
            ! j0 = j0 + 1
          ! end do
         end if
         return
         end if  ! end of if flat barrier case
   !=========================================================
   !                     |
   !                     |  vertical case
   !                     |
   !=========================================================
        if (vertical.eqv..true.) then
          if (i_0 /= i_e) then
            write(*,*) "The wall is not vertical"
            return
          end if

           ! distance away from physical edge
           alpha = abs(x_0 - xe(i_0+1))/dx
           if (alpha.lt.1d-7) then
             alpha = 0.d0
           end if

           ! the right h-boxes
           if (i_0 == size(xe)) then
             write(*,*) "Wall is too left, h-box will be off the grid"
             return
           end if
           if (i_0 == size(xe)-1) then
             write(*,*) "Can only make one layer of h-box"
             if (layer==2) then
               return
             end if
           end if
           i_2 = i_0 + 2
           i_1 = i_0 + 1
           i_3 = i_0 + 3
           if (layer == 1) then
             j0 = j_0
             ! do j = 1,m
       q_hbox(1:3,1-mbc:mx+mbc) = (1-alpha)*qold(1:3,i_1,:)+alpha*qold(1:3,i_2,:)
       aux_hbox(1,1-mbc:mx+mbc) = (1-alpha)*aux(1,i_1,:)+alpha*aux(1,i_2,:)
               ! j0 = j0 + 1
             ! end do
           end if
           if (layer == 2) then
             j0 = j_0
             ! do j = 1,m
       q_hbox(1:3,1-mbc:mx+mbc) = (1-alpha)*qold(1:3,i_2,:)+alpha*qold(1:3,i_3,:)
       aux_hbox(1,1-mbc:mx+mbc) = (1-alpha)*aux(1,i_2,:)+alpha*aux(1,i_3,:)
               ! j0 = j0 + 1
             ! end do
           end if
           return
         end if ! end of vertical wall case
       end subroutine up_hboxes

end module hbox2D

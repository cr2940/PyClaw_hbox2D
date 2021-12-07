! CODES FOR SRD:


module aux_module_SRD2

  ! use hbox_layer
  use redistribute2D
  implicit none
  real(8), parameter:: tol = 5d-13, pi=3.14159265359d0

contains
  function dist(A,B) result(distance)
    implicit none
    real(kind=8) :: A(:),B(:)
    integer :: n ,i
    real(kind=8) :: distance
    n = size(A)
    distance = 0.d0
    do i=1,n
      distance = distance + (A(i)-B(i))**2
    end do
    distance = sqrt(distance)
    ! distance = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
  end function dist
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

  subroutine hbox_SRD_update(mx,q_hbox,area_hboxes,N_hbox,whose)
    ! updates hboxes with the final step of SRD'ing in the transverse direction
    implicit none
    real(8) :: q_hbox(:,:), area_hboxes(:)
    integer :: N_hbox(:), whose(:) ! whose is the index pointer of neighbor for hboxes.
    integer :: i,N,mx
    real(8) :: nhood_vals(3,mx*4)
    N = size(q_hbox,2)
    ! print *, "N:",N
    do i = 1,N
      ! print *, "whose:",whose(i)
      if (whose(i).eq.i) then
        nhood_vals(:,i) = q_hbox(:,i)
        ! print *, "Nhood: itself", nhood_vals(:,i)
      else
        nhood_vals(:,i) = (area_hboxes(i)/N_hbox(i)+area_hboxes(whose(i))/N_hbox(whose(i)))**(-1)&
        * (area_hboxes(i)/N_hbox(i)*q_hbox(:,i)+area_hboxes(whose(i))/N_hbox(whose(i))*q_hbox(:,whose(i)))
        ! print *, "Nhood: not itself", nhood_vals(:,i)

      end if
    end do
    do i=1,N
      if (N_hbox(i).eq.1) then ! the case it is not overlapped by anyone
        q_hbox(:,i) = nhood_vals(:,i)
        ! print *, "result itself:", q_hbox(:,i)
      else
        q_hbox(:,i) = nhood_vals(:,i)/N_hbox(i) + nhood_vals(:,i-1)/N_hbox(i)
        ! print *, "result not itself:", q_hbox(:,i)
      end if
    end do
  end subroutine

  subroutine conserve_correct(q_hbox,num_frags,index_frags,area_frags,cov_inds,qinc,hbox_area,mx,my,mbc,dx,dy,areaij)
    ! corrects for conservation of mass with the hbox method by going through the affected indices,
    ! and if not fully covered, takes portion of its increment and subtract
    implicit none
    real(8) :: q_hbox(:,:),area_frags(:,:),dx,dy,hbox_area(:)
    real(8) :: areaij(:,:)
    integer :: num_frags(:), index_frags(:,:,:), cov_inds(:,:)
    integer :: m,i,j,look(2),mx,mbc,my
    real(8) :: qinc(3,1-mbc:mx+mbc,1-mbc:my+mbc)
    m = size(q_hbox,2)
    do i = 1,m
      do j=1,num_frags(i)
        look = index_frags(:,j,i)
        ! print*, qinc
        if (.not. ind_in_inds(look,cov_inds)) then
          q_hbox(:,i) = q_hbox(:,i) + qinc(:,look(1),look(2))*area_frags(j,i)*hbox_area(i)/&
              (areaij(look(1),look(2)))

        end if
        ! aux_hbox_d(i) = aux_hbox_d(i) + aux(1,look(1),look(2))*area_frags_d(j,i)
      end do
    enddo
  end subroutine
  subroutine rotate_qhbox_orig(N,q_hbox,q_hbox_r,intersections,ot)
    implicit none
    real(8) :: q_hbox(:,:), intersections(:,:)
    logical :: ot
    real(8) :: x,y,n_vec(2),n_vec2(2),t_vec(2),t_vec2(2)
    real(8) :: vec1(2),vec2(2),q_hbox_r(3,N)
    integer :: i , N
    N = size(q_hbox,2)
    do i=1,N
      x = intersections(1,i+1)- intersections(1,i)
      y = intersections(2,i+1) - intersections(2,i)
      n_vec = -(/-y,x/)
      n_vec2 = (/y,x/) ! for turning back to original
      n_vec = n_vec/(sqrt(x**2+y**2))
      n_vec2 = n_vec2/(sqrt(x**2+y**2))
      t_vec = -(/x,y/)
      t_vec2 = (/-x,y/) ! for turning back to orignal
      t_vec = t_vec/(sqrt(x**2+y**2))
      t_vec2 = t_vec2/sqrt(x**2+y**2)
      if (.not. ot) then
        vec1 = t_vec
        vec2 = -n_vec
      else
        vec1 = n_vec2
        vec2 = t_vec2
      end if
      call rotate_state(q_hbox(:,i),q_hbox_r(:,i),vec1,vec2)
    end do
  end subroutine

  subroutine grid_update_hbox(mx,my,mbc,dx,dy,q_hbox_down,q_hbox_up,qnew,qnew2,uniq_indices_down,uniq_indices_up,&
    indices_data_area_down,indices_data_area_up,un_area_ij,up_area_ij,count_down,count_up,&
    cov_ind_down,cov_ind_up,hbox_area_down,hbox_area_up,indices_data_up,indices_data_down)
    ! here the q_hbox is all updated ready to be applied to the underlying grid
    ! also qnew, qnew2 are edited using small cell update subroutine
    ! so this subroutine takes the updated hbox values and applies them to the grid
    ! loop through the affected cells
    ! look at which hbox covers it and how much
    ! and if fills the whole cell, then done, but if not, take the remaining area and weight the cut cell godunov update and average that in the new update
    implicit none
    integer ::  mx,my,mbc
    real(8) :: q_hbox_down(:,:), q_hbox_up(:,:),dx,dy
    real(8) :: un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(8) :: qnew(3,1-mbc:mx+mbc,1-mbc:my+mbc), qnew2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: uniq_indices_down(:,:), uniq_indices_up(:,:), count_down(:),count_up(:)
    real(8) :: indices_data_area_down(:,:), indices_data_area_up(:,:), q_orig(3)
    real(8) :: area_sum, hbox_area_down(:), hbox_area_up(:), area_actual,rem,q_hb(3)
    integer :: i,n1,n2, ind(2),j,k,cov_ind_down(:,:),cov_ind_up(:,:)
    integer :: indices_data_up(:,:), indices_data_down(:,:)

    ! for upper grid
    n1 = size(uniq_indices_up,2)
    do i =1,n1
      ind = uniq_indices_up(:,i)
      area_actual = up_area_ij(ind(1),ind(2))*dx*dy
      q_orig = qnew2(:,ind(1),ind(2))
      ! print*, "ind:",ind
      ! print*, "j :", count_up(i)
      ! print*, "qn2:", qnew2(:,ind(1),ind(2))
      j = count_up(i)
      if (ind_in_inds(ind,cov_ind_up)) then ! the affected cell fully covered by hboxes case
        j = count_up(i)
        qnew2(:,ind(1),ind(2)) = 0.d0
        do k=1,j
           qnew2(:,ind(1),ind(2)) =qnew2(:,ind(1),ind(2)) + hbox_area_up(indices_data_up(k,i))*&
             indices_data_area_up(k,i)/area_actual * q_hbox_up(:,indices_data_up(k,i)) !
        end do

      else ! the affected cell partially covered by hboxes case
         ! print *, "Qorig:", q_orig
         q_hb = 0.d0
         area_sum=0.d0
         do k=1,j
           ! if (indices_data_up(k,i).eq.0) then
           !   cycle
           ! end if
           ! print *, "Inddata:",indices_data_up(k,i)
           area_sum =area_sum + hbox_area_up(indices_data_up(k,i))*&
             indices_data_area_up(k,i)!/area_actual
             ! if (area_sum .gt. area_actual) then
             !   area_sum =area_actual !area_sum - hbox_area_up(indices_data_up(k,i))*&
             !     ! indices_data_area_up(k,i)
             ! end if
             ! print *, "coeff",hbox_area_up(indices_data_up(k,i))*indices_data_area_up(k,i)/area_actual
             ! print *, "qbox",q_hbox_up(:,indices_data_up(k,i))

           q_hb = q_hb + hbox_area_up(indices_data_up(k,i))*&
             indices_data_area_up(k,i)/area_actual* q_hbox_up(:,indices_data_up(k,i)) !
         end do
         rem = (area_actual - area_sum)/area_actual
         ! print *, "areas: u",area_actual, area_sum, rem
         qnew2(:,ind(1),ind(2)) = max(rem,0.d0)*q_orig + q_hb !/area_actual
       end if
       if (dist(q_orig,qnew2(:,ind(1),ind(2))) .lt. 1d-2) then
         qnew2(:,ind(1),ind(2)) = q_orig
       end if
       ! print*, "qn2 after: ", qnew2(:,ind(1),ind(2))
       ! print *, "qhbox:",q_hb,rem

       ! print*, "qn2 after:", qnew2(:,ind(1),ind(2))!, ind_in_inds(ind,cov_ind_up)

     end do
     ! for lower grid
     n2 = size(uniq_indices_down,2)
     do i =1,n2
       ind = uniq_indices_down(:,i)
       area_actual = un_area_ij(ind(1),ind(2))*dx*dy
       q_orig = qnew(:,ind(1),ind(2))
       j = count_down(i)
       ! print*, "qn:", qnew(:,ind(1),ind(2))

       if (ind_in_inds(ind,cov_ind_down)) then ! the affected cell fully covered by hboxes case
         j = count_down(i)
         qnew(:,ind(1),ind(2)) = 0.d0
         do k=1,j
            qnew(:,ind(1),ind(2)) =qnew(:,ind(1),ind(2)) + hbox_area_down(indices_data_down(k,i))*&
              indices_data_area_down(k,i)/area_actual * q_hbox_down(:,indices_data_down(k,i)) !
         end do

       else ! the affected cell partially covered by hboxes case
          ! print *, "Qorig:", q_orig
          q_hb = 0.d0
          area_sum=0.d0
          do k=1,j
            if (indices_data_down(k,i).eq.0) then
              cycle
            end if
            area_sum =area_sum + hbox_area_down(indices_data_down(k,i))*&
              indices_data_area_down(k,i)!/area_actual
              ! if (area_sum .gt. area_actual) then
              !   area_sum =area_actual !area_sum - hbox_area_up(indices_data_up(k,i))*&
              !     ! indices_data_area_up(k,i)
              ! end if
            q_hb = q_hb + hbox_area_down(indices_data_down(k,i))*&
              indices_data_area_down(k,i)/area_actual * q_hbox_down(:,indices_data_down(k,i))!
          end do
          rem = (area_actual - area_sum)/area_actual
          ! print *, "areas: d",area_actual, area_sum, rem
            qnew(:,ind(1),ind(2)) = max(rem,0.d0)*q_orig + q_hb !/area_actual
        end if
        if (dist(q_orig,qnew(:,ind(1),ind(2))) .lt. 1d-2) then
          qnew(:,ind(1),ind(2)) = q_orig
        end if
        ! print*, "qn after:", qnew(:,ind(1),ind(2))!, ind_in_inds(ind,cov_ind_up)

        ! print *, "qhbox:",q_hb,rem
        ! print*, "qn after:", qnew(:,ind(1),ind(2))!, ind_in_inds(ind,cov_ind_up)

      end do
    end subroutine


    ! determine neighborhood info for the hboxes.
    !! UPDATE::

    ! for each AFFECTED FRAGMENT INDICES , KEEP
    ! copy of indices of HBOXES that cover it, the AREAS of those hboxes,
    ! TALLY of how many indices, If they dont add to the whole area of cell, take remainder and use original godunov update. (you have to use small cells update)

    ! q_hbox(:,i) = q_hbox(:,i) - dtdx*(amdq_wall)
    ! call rn2(q_hbox,q_hbox) --> apdq,amdq
    ! q_hbox(:,i) = q_hbox(:,i)- dtdx/area(i) * (apdq(:,i) +amdq(:,i))
    !
    ! go through each hbox, get/determine N_ij
    ! go through each hbox, point to the cell that is its neighboring component
    ! go through each hbox, get the covering neighborhood value
    ! go through each hbox, divide the neighborhood value by the N_ij


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


    subroutine small_cells_update_correct(ii,jj,dtdx,qold,qold2,qnew,qnew2,aux,&              !  [   ]
      intersections,type_supper,type_sunder,lengths_supper,lengths_sunder,&
        area_sunder,area_supper,ixy,N_cells,xlower,ylower,xupper,yupper,&
        mx,my,mbc,x_0,y_0,x_e,y_e,wall_height,fm,fp,gm,gp,fm2,fp2,gm2,gp2,ot,&
        q_hbox_d,q_hbox_u,hbox_areas_d,hbox_areas_u,comp_cov_inds_d,comp_cov_inds_u)
      ! takes qnew/qnew2 and uses qold/qold2 to correct the update so that its specific to the small cells' geometry
      ! ie computes Q^ in the paper for both up and under small cells
      implicit none

      ! input / output
      integer :: ii(:),jj(:),type_supper(:),type_sunder(:),ixy,N_cells,mx,my,mbc
      integer :: comp_cov_inds_d(:,:),comp_cov_inds_u(:,:)
      real(8):: qold(3,1-mbc:mx+mbc,1-mbc:my+mbc),qold2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) :: qnew(3,1-mbc:mx+mbc,1-mbc:my+mbc),qnew2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) :: aux(1,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) :: intersections(:,:),lengths_sunder(:,:),lengths_supper(:,:)
      real(8) :: area_supper(:),area_sunder(:),xlower,ylower,xupper,yupper
      real(8) :: x_0,y_0,x_e,y_e,wall_height,dtdx
      real(8) fp(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) fm(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) gp(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) gm(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) fp2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) fm2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) gp2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(8) gm2(3,1-mbc:mx+mbc,1-mbc:my+mbc)

      ! local
      real(8) :: n_vec(2),t_vec(2),x,y,n_vec2(2),t_vec2(2),vec1(2),vec2(2)
      real(8) :: q_hbox_d(3,mx*4), q_hbox_u(3,mx*4), aux_hbox_u(mx*4),aux_hbox_d(mx*4)
      logical :: L2R, R2L
      integer :: is,js,i,j,k,m,num_frags_d(mx*4),index_frags_d(2,4,mx*4)
      integer :: num_frags_u(mx*4),index_frags_u(2,4,mx*4),look(2),ll,mm
      real(8) :: hbox_areas_u(mx*4), area_frags_u(4,mx*4)
      real(8) :: hbox_areas_d(mx*4), area_frags_d(4,mx*4)
      real(8) :: wave_wall(3,3), amdq_wall(3),apdq_wall(3),s_wall(3)
      real(8) :: hL,hR,huL,huR,hvL,hvR,bL,bR,hstarL,hstarR,ustarL,ustarR
      logical :: lexist,ot
      real(8) :: dir1(3),dir2(3),dir3(3),dir4(3),s(3),fwave(3,3),apdq(3)
      real(8) :: qtemp1(3),qtemp2(3),qtemp3(3)

      ! get the hboxes for normal barrier flux calculation:
      inquire (file="./hbox_data2.txt",exist=lexist)
      ! print *, "EXIST:",lexist
      open (unit=2,file="./hbox_data2.txt",FORM="FORMATTED",STATUS="OLD",&
      ACTION="READ",access='sequential')
      rewind 2
      read(2,*) m
      read(2,*) num_frags_d
      read(2,*) num_frags_u
      read(2,*) index_frags_d
      read(2,*) index_frags_u
      read(2,*) hbox_areas_d
      read(2,*) hbox_areas_u
      read(2,*) area_frags_u
      read(2,*) area_frags_d
      close(2)!,status="keep")
  100 continue
      q_hbox_d = 0.d0
      q_hbox_u = 0.d0
      ! q_hbox_d(1,:) = 1.d0
      ! q_hbox_u(1,:) = 1.d0
      aux_hbox_d = -2.d0
      aux_hbox_u = -2.d0

      do i = 1,m
        do j=1,num_frags_d(i)
          look = index_frags_d(:,j,i)
          ! print *, "THE VLAUE:", qold(:,look(1),look(2))
          q_hbox_d(:,i) = q_hbox_d(:,i) + qold(:,look(1),look(2))*area_frags_d(j,i)
          ! aux_hbox_d(i) = aux_hbox_d(i) + aux(1,look(1),look(2))*area_frags_d(j,i)
        end do
      enddo
      do i = 1,m
        do j=1,num_frags_u(i)
          look = index_frags_u(:,j,i)
          ! print *, "THE VLAUE 2:", qold2(:,look(1),look(2))
          q_hbox_u(:,i) = q_hbox_u(:,i) + qold2(:,look(1),look(2))*area_frags_u(j,i)
          ! aux_hbox_u(i) = aux_hbox_u(i) + aux(1,look(1),look(2))*area_frags_u(j,i)
        end do
      enddo

       ot = .false.
      ll = 1  ! lower
      mm = 1 ! upper
      ! if (ixy.eq.1) then

        do i = 1,N_cells
          is = ii(i)
          js = jj(i)
          x = intersections(1,i+1)- intersections(1,i)
          y = intersections(2,i+1) - intersections(2,i)
          ! print *, "X: ", intersections(1,i), intersections(1,i+1)
          ! print *, "Y: ", intersections(2,i), intersections(2,i+1)
          ! print *, "i:  ", i , "slope: ", y/x
          n_vec = -(/-y,x/)
          n_vec2 = (/y,x/) ! for turning back to original
          n_vec = n_vec/(sqrt(x**2+y**2))
          n_vec2 = n_vec2/(sqrt(x**2+y**2))
          t_vec = -(/x,y/)
          t_vec2 = (/-x,y/) ! for turning back to orignal
          t_vec = t_vec/(sqrt(x**2+y**2))
          t_vec2 = t_vec2/sqrt(x**2+y**2)

          call rotate_state(q_hbox_d(:,i),qtemp1, &
          -ll*n_vec,-ll*t_vec)
          call rotate_state(q_hbox_u(:,i),qtemp2,&
           -mm*n_vec,-mm*t_vec)

          hL = qtemp2(1)
          hR = qtemp1(1)
          huL = qtemp2(2)
          huR = qtemp1(2)
          hvL= qtemp2(3)
          hvR = qtemp1(3)
          bL = aux_hbox_u(i)
          bR = aux_hbox_d(i)

          call barrier_passing(hL,hR,huL,huR,bL,bR,wall_height,&
                    L2R,R2L,hstarL,hstarR,ustarL,ustarR)
          if (L2R .or. R2L) then
            ot = .true.
          end if
          call redistribute_fwave(1,qtemp2,qtemp1,aux_hbox_u(i),&
             aux_hbox_d(i),wall_height,1,wave_wall,s_wall,amdq_wall,apdq_wall,3,&
             3,L2R,R2L)

          call rotate_state(amdq_wall,amdq_wall,-t_vec,-n_vec2)
          call rotate_state(apdq_wall,apdq_wall,-t_vec,-n_vec2)

         qtempu(:,i) = q_hbox_u(:,i) - dtdx*amdq_wall
         qtempd(:,i) = q_hbox_d(:,i) - dtdx*apdq_wall


        if (.not. ind_in_inds((/is,js/),comp_cov_inds_d)) then
          select case (type_sunder(i))
          case (1) ! TYPE 1 Cut cell
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(2,i)*&
            (fm(:,is+1,js))!+f(qold(:,is,js),1))
            !call rotate_state(qold(:,is,js),qold(:,is,js),ll*n_vec,ll*t_vec)
            ! print *, "ROT LOW:", qold(:,is,js)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              ( apdq_wall )!f(qold(:,is,js),1)
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            qnew(:,is,js) = qnew(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(3,i)*(gp(:,is,js))!&

            ! print *, "QNEW:", qnew(:,is,js)

          case (2) ! TYPE 2 Cut cell
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*fm(:,is+1,js)
            !call rotate_state(qold(:,is,js),qold(:,is,js),mm*n_vec,mm*t_vec)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              (apdq_wall) !f(qold(:,is,js),1)
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            qnew(:,is,js) = qnew(:,is,js) -dtdx/area_sunder(i)*((lengths_sunder(2,i))*gm(:,is,js+1)+&
               (lengths_sunder(4,i)) * gp(:,is,js))! +dtdx*(lengths_sunder(2,i)-lengths_sunder(4,i))*&
                 ! f(qold(:,is,js),2)
          case (3,6) ! TYPE 3/6 Cut cell
            ! under small cell
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*((lengths_sunder(2,i))*fm(:,is+1,js) + &
              (lengths_sunder(4,i))*fp(:,is,js)) !- dtdx/area_sunder(i)*(lengths_sunder(2,i)-lengths_sunder(4,i))&
               ! * (f(qold(:,is,js),1))
           qnew(:,is,js) = qnew(:,is,js) - dtdx/area_sunder(i)*(gp(:,is,js)) !-f(qold(:,is,js),2)+
           !call rotate_state(qnew(:,is,js),qnew(:,is,js),ll*n_vec,ll*t_vec)
           qnew(:,is,js) = qnew(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
             (apdq_wall)!f(qold(:,is,js),1)
           !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)

          case (4) ! TYPE 4 Cut cell
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*(lengths_sunder(5,i))*(fp(:,is,js))&!-f(qold(:,is,js),1))
             - dtdx/area_sunder(i)*(fm(:,is+1,js))!+f(qold(:,is,js),1))
            !call rotate_state(qold(:,is,js),qold(:,is,js),ll*n_vec,ll*t_vec)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              (apdq_wall)!f(qold(:,is,js),1)
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            qnew(:,is,js) = qnew(:,is,js) - dtdx/area_sunder(i)*(lengths_sunder(2,i))*& !f(qold(:,is,js),2) +
              (gm(:,is,js+1)) - dtdx/area_sunder(i)*(gp(:,is,js))!-f(qold(:,is,js),2))!

          case (5) ! TYPE 5 Cut cell
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*(lengths_sunder(2,i))*(fm(:,is+1,js))&
               -dtdx/area_sunder(i)*fp(:,is,js)
            !call rotate_state(qold(:,is,js),qold(:,is,js),ll*n_vec,ll*t_vec)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              (apdq_wall)!f(qold(:,is,js),1) -
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            ! f(qold(:,is,  js),1) +
            qnew(:,is,js) = qnew(:,is,js) -dtdx/area_sunder(i)*(lengths_sunder(5,i))*&
              (gm(:,is,js+1)) -dtdx/area_sunder(i)*gp(:,is,js)  !f(qold(:,is,js),2) +
          case (7)
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*fp(:,is,js)
            !call rotate_state(qold(:,is,js),qold(:,is,js),ll*n_vec,ll*t_vec)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              (apdq_wall)!f(qold(:,is,js),1)
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            qnew(:,is,js) = qnew(:,is,js)- dtdx/area_sunder(i)*((lengths_sunder(4,i))*gm(:,is,js+1) &
              + (lengths_sunder(2,i))*gp(:,is,js)) !+dtdx*(lengths_sunder(4,i)-lengths_sunder(2,i))&
                ! * (f(qold(:,is,js),2))
          case(8)
            ! under small cell
            qold(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(3,i)*(fp(:,is,js))
            !call rotate_state(qold(:,is,js),qold(:,is,js),ll*n_vec,ll*t_vec)
            qnew(:,is,js) = qold(:,is,js) - dtdx/area_sunder(i)*lengths_sunder(1,i)*&
              ( apdq_wall)!f(qold(:,is,js),1)
            !call rotate_state(qnew(:,is,js),qnew(:,is,js),vec1,vec2)
            ! f(qold(:,is,js),1)&
            qnew(:,is,js) = qnew(:,is,js) - dtdx/area_sunder(i)*(lengths_sunder(2,i))*&
              (gp(:,is,js)) !f(qold(:,is,js),2)

          end select
      end if

       if (.not. ind_in_inds((/is,js/),comp_cov_inds_u)) then
          select case (type_supper(i))
          case (1) ! TYPE 1 Cut cell
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*(lengths_supper(2,i))*&
             (fm2(:,is+1,js)) - dtdx/area_supper(i)*(fp2(:,is,js))! & + f(qold2(:,is,js),1)
              ! - f(qold2(:,is,js),1))
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            ! print *, "ROT: ", qold2(:,is,js)
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall)!f(qold2(:,is,js),1)+
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            ! print *, "QNEW2:", qnew2(:,is,js)
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*(lengths_supper(5,i))&
            *(gp2(:,is,js) ) -dtdx/area_supper(i)*(gm2(:,is,js+1))!& !- f(qold2(:,is,js),2)
            !f(qold2(:,is,js),2))!

          case (2) ! TYPE 2 Cut cell
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*fp2(:,is,js)
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              ( amdq_wall)!f(qold2(:,is,js),1)+
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            qnew2(:,is,js) = qnew2(:,is,js) -dtdx/area_supper(i)*((lengths_supper(2,i))*gm2(:,is,js+1)+&
               (lengths_supper(4,i)) * gp2(:,is,js)) !+dtdx*(lengths_supper(2,i)-lengths_supper(4,i))*&
                 ! f(qold2(:,is,js),2)

          case (3,6)
            ! upper small cell
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*((lengths_supper(2,i))*fm2(:,is+1,js) + &
              (lengths_supper(4,i))*fp2(:,is,js)) !- dtdx/area_supper(i)*(lengths_supper(2,i)-lengths_supper(4,i))&
                 ! * (f(qold2(:,is,js),1))
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*(gm2(:,is,js+1))!f(qold2(:,is,js),2)+
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall) !f(qold2(:,is,js),1) +
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)

          case (4)
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(3,i)*(fp2(:,is,js))!&
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall)!f(qold2(:,is,js),1) +
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            !! -f(qold2(:,is,js),1))
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*lengths_supper(2,i)* &
              (gm2(:,is,js+1))!f(qold2(:,is,js),2) +

          case (5)
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(2,i)*(fm2(:,is+1,js))
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall)!f(qold2(:,is,js),1) +
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            ! f(qold2(:,is,js),1) +
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*lengths_supper(3,i)* &
              (gm2(:,is,js+1))!f(qold2(:,is,js),2) +

          case (7)
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*fp2(:,is,js)
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall) !f(qold2(:,is,js),1) +
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            qnew2(:,is,js) = qnew2(:,is,js)- dtdx/area_supper(i)* ((lengths_supper(4,i))*gm2(:,is,js+1) &
              + (lengths_supper(2,i))*gp2(:,is,js))!+dtdx*(lengths_supper(4,i)-lengths_supper(2,i))&
                ! * (f(qold2(:,is,js),2))

          case(8)
            ! upper small cell
            qold2(:,is,js) = qold2(:,is,js) - dtdx/area_supper(i)*(lengths_supper(5,i))*(fp2(:,is,js))&
             - dtdx/area_supper(i)*fm2(:,is+1,js)
            !call rotate_state(qold2(:,is,js),qold2(:,is,js),mm*n_vec,mm*t_vec)
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*lengths_supper(1,i)*&
              (amdq_wall)!f(qold2(:,is,js),1) +
            !call rotate_state(qnew2(:,is,js),qnew2(:,is,js),vec1,vec2)
            ! f(qold2(:,is,&  js),1)
            qnew2(:,is,js) = qnew2(:,is,js) - dtdx/area_supper(i)*(lengths_supper(2,i))*&
              (gp2(:,is,js)) - dtdx/area_supper(i)*gm2(:,is,js+1) !f(qold2(:,is,js),2)

          end select
         end if
          ! ! print*, "*******TYPE: ", type_supper(i), "***********"
          ! !
          ! ! print*, "the undersmall cell: post" , qnew(:,is,js)
          ! ! print*, "the uppersmall cell: post", qnew2(:,is,js)
        end do
       

    end subroutine


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
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.tol))then
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
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.tol))then
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
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.tol))then
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
        if (area_sunder(i) .lt. 0.5d0 .and. abs(area_sunder(i)-0.5d0).gt.tol) then
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
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.tol))then
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
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.tol))then
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
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.tol))then
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
        if (area_supper(i) .lt. 0.5d0.and. abs(area_supper(i)-0.5d0).gt.tol) then
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


    subroutine area_cells(ii,jj,mx,my,mbc,area_sunder,area_supper,up_area_ij,un_area_ij)
      implicit none
      integer :: ii(:),jj(:),mx,my,mbc
      real(8) :: area_sunder(:),area_supper(:)
      real(8),intent(out) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
      ! local
      integer :: N,i

      ! initialization of area matrix
      up_area_ij = 1.d0
      un_area_ij = 1.d0

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
      real(8) :: beta,alpha,nhood_val(3),val(3)
      integer :: num_ij2

      ! DO UPPER SIDE FIRST:
      N = size(ii)! size(all_undercells_i)
      k=1
      do i=1,N
        i0 = ii(i)!all_uppercells_i(i)
        j0 = jj(i)!all_uppercells_j(i)
        ! print *, "Q UP BEFORE: ", qnew2(:,i0,j0)

        ! if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
        !   array2(i) = k
        !   k =k +1
        ! else
        !   array2(i)=k-1
        ! end if
        ! print*, "ARRAY2:", array2(i)
        num_ij = N_ij_up(i0,j0)
        alpha = up_area_ij(i0,j0)
        ! print*, "AREA: ",alpha
        ! print *, "NUM_IJ:", num_ij
        ! if (num_ij .eq. 2) then
        !   qnew2(:,i0,j0) = qnew2(:,i0,j0)/2.d0
        if (alpha.lt.0.5d0) then
          ! print *, "area: UP", up_area_ij(i0,j0)
          ! get neighborhood value:
          ! if (up_area_ij(i0,j0).lt.0.5d0) then!  .and. abs(up_area_ij(i0,j0)-0.5d0).gt.1d-10)then
            i_n = i0!upS_cells_i(array2(i))
            j_n = j0+1!upS_cells_j(array2(i))
            ! print *, "In,jn:",i_n,j_n
            ! print *, "Q UP NEIGHBOR BEFRE: ", qnew2(:,i_n,j_n)
            ! print *, "BETA: " , beta
            ! print *, "ALPAH: ", alpha
            ! print *, "CELL INTERST VAL: ", qnew2(:,i0,j0)
            ! print *, "NEIGHOR CELL VAL: ", qnew(:,i_n,j_n)
            beta = up_area_ij(i_n,j_n)
            ! print*, "Neighbor AREA: ",beta
            ! print *, "N_ij:", N_ij_up(i_n,j_n)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew2(:,i0,j0)+&
                 beta/2.d0*qnew2(:,i_n,j_n))
            ! print *, 'Nhood val:', nhood_val(1)
            qnew2(:,i0,j0) = nhood_val
            num_ij2 = N_ij_up(i_n,j_n)
            ! if (num_ij2 .eq. 2) then
              val =qnew2(:,i_n,j_n)/2.d0
              ! print *, "val'", val
              qnew2(:,i_n,j_n) = val + nhood_val/2.d0
              ! print *, "averag val:", qnew2(1,i_n,j_n)
            ! end if
          end if
        ! end if
        ! print *, "Q UP AFTER: ", qnew2(:,i0,j0)
        ! print *, "Q  UP NEIGHBOR AFTER: ", qnew2(:,i_n,j_n)
      end do

      ! DO UNDER SIDE NEXT:
      N = size(ii)! size(all_undercells_i)
      k = 1
      do i=1,N
        ! print *, "Q BEFORE: ", qnew(:,i0,j0)
        i0 = ii(i)!all_undercells_i(i)
        j0 = jj(i)!all_undercells_j(i)
        ! if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
        !   array(i) = k
        !   k = k +1
        ! else
        !   array(i)=k-1
        ! end if
        ! print *, " I0,J0:", i0,j0
        ! print *, "ARRAY:", array(i)
        num_ij = N_ij_un(i0,j0)
        alpha = un_area_ij(i0,j0)
        ! print*, "AREA: ",alpha
        ! print *, "NUM_IJ:", num_ij
        ! if (num_ij .eq. 2) then
        !   qnew(:,i0,j0) = qnew(:,i0,j0)/2.d0
        if (alpha.lt.0.5d0) then
          ! get neighborhood value:
          ! print *, "area: UN", un_area_ij(i0,j0)
          ! if (un_area_ij(i0,j0).lt.0.5d0) then! .and. abs(un_area_ij(i0,j0)-0.5d0).gt.1d-10)then
            i_n = i0!unS_cells_i(array(i))
            j_n = j0-1!unS_cells_j(array(i))
            ! print*, "IN,JN:", i_n,j_n
            ! print *, "Q NEIGHBOR BEFRE: ", qnew(:,i_n,j_n)
            beta = un_area_ij(i_n,j_n)
            ! print*, "Neighbor AREA: ",beta
            ! print *, "N_ij:", N_ij_up(i_n,j_n)
            ! alpha = un_area_ij(i0,j0)
            ! print *, "BETA: " , beta
            ! print *, "ALPAH: ", alpha
            ! print *, "CELL INTERST VAL: ", qnew2(:,i0,j0)
            ! print *, "NEIGHOR CELL VAL: ", qnew(:,i_n,j_n)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew(:,i0,j0)+&
                 beta/2.d0*qnew(:,i_n,j_n))
           ! print *, 'Nhood val:', nhood_val(1)
            qnew(:,i0,j0) = nhood_val
            ! print *, "NHOOD:", nhood_val
            num_ij2 = N_ij_un(i_n,j_n)
            ! if (num_ij2 .eq. 2) then
              val =qnew(:,i_n,j_n)/2.d0
              ! print *, "val'", val
              qnew(:,i_n,j_n) = val + nhood_val/2.d0
              ! print *, "averag val:", qnew2(1,i_n,j_n)
            ! end if
          end if
        ! end if
        ! print *, "Q AFTER: ", qnew(:,i0,j0)
        ! print *, "Q NEIGHBOR AFTER: ", qnew(:,i_n,j_n)
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

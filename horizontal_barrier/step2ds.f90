!     ==========================================================
    subroutine step2ds(maxm,num_eqn,num_waves,num_aux,num_ghost,mx,my, &
                        qold,qnew,aux,dx,dy,dt,method,mthlim,cfl, &
                        qadd,fadd,gadd,q1d,dtdx1d,dtdy1d, &
                        aux1,aux2,aux3,work,mwork,ids,use_fwave,rpn2,rpt2,vert,horz,&
                        bar_index_i,bar_ht,bar_loc)
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
    use hbox2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit double precision (a-h,o-z)
    double precision :: qold(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qnew(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision ::  q1d(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: qadd(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: fadd(num_eqn, 1-num_ghost:maxm+num_ghost)
    double precision :: gadd(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
    double precision :: aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: aux1(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux2(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux3(num_aux, 1-num_ghost:maxm+num_ghost)

    double precision :: dtdx1d(1-num_ghost:maxm+num_ghost)
    double precision :: dtdy1d(1-num_ghost:maxm+num_ghost)
    integer :: method(7),mthlim(num_waves)
    logical ::          use_fwave
    double precision :: work(mwork)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dimension amdqd2(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension waved2(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    dimension sd2(num_waves,1-num_ghost:maxm+num_ghost)
    dimension apdqd2(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension amdqd1(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension waved1(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    dimension sd1(num_waves,1-num_ghost:maxm+num_ghost)
    dimension apdqd1(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension amdqu1(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension waveu1(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    dimension su1(num_waves,1-num_ghost:maxm+num_ghost)
    dimension apdqu1(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension amdqu2(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension waveu2(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    dimension su2(num_waves,1-num_ghost:maxm+num_ghost)
    dimension apdqu2(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension apdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension amdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    dimension s_r(3),fwave_r(3,3)
    double precision :: amdq_hbox(num_eqn,4),apdq_hbox(num_eqn,4)
    real(kind=8):: q_hbox_d1(3,1-num_ghost:maxm+num_ghost),q_hbox_d2(3,1-num_ghost:maxm+num_ghost)
    real(kind=8):: aux_hbox_d1(1,1-num_ghost:maxm+num_ghost),q_hbox_u1(3,1-num_ghost:maxm+num_ghost)
    real(kind=8):: aux_hbox_d2(1,1-num_ghost:maxm+num_ghost),aux_hbox_u1(1,1-num_ghost:maxm+num_ghost)
    real(kind=8) :: aux_hbox_u2(1,1-num_ghost:maxm+num_ghost)
    real(kind=8) :: q_hbox_u2(3,1-num_ghost:maxm+num_ghost)
    ! barrier parameters
    real(kind=8) :: bar_l
    real(kind=8)::x_0,x_e,y_0,y_e,alpha,bar_height,bar_ht,bar_loc
    integer :: i_0,i_e,j_0,j_e,bar_index,bar_index_i
    logical :: horz, vert
    double precision :: temp(3,1-num_ghost:mx+num_ghost)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external :: rpn2,rpt2
    ! barrier param set on setup python file
    ! common /cparam/ bar_index_i,bar_ht,bar_loc
!f2py intent(out) cfl
!f2py intent(in,out) qnew
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


    index_capa = method(6)
    num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! hboxes:
    bar_index= bar_index_i
    bar_height = bar_ht
    bar_l = bar_loc
    ! write(*,*) "Y0:", bar_l,bar_index,bar_height
    if (horz) then
      y_0 = bar_l
      x_0 = 0.d0 - 2*dx
      x_e = 1.d0 + 2*dx
      j_0 = bar_index
      i_0 = 1-num_ghost
      y_e = bar_l
      j_e = bar_index
      i_e = mx+num_ghost
    else if (vert) then
      x_0 = bar_l
      y_0 = 0.d0 - 2*dx
      y_e = 1.d0 + 2*dx
      i_0 = bar_index
      j_0 = 1-num_ghost
      x_e = bar_l
      i_e = bar_index
      j_e = mx+num_ghost
    end if
    call down_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,2,&
        qold,aux,num_eqn,1,num_ghost,mx,my,maxm,dx,dy,q_hbox_d2,aux_hbox_d2,alpha)
    call down_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,1,&
        qold,aux,num_eqn,1,num_ghost,mx,my,maxm,dx,dy,q_hbox_d1,aux_hbox_d1,alpha)
    call up_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,1,&
        qold,aux,num_eqn,1,num_ghost,mx,my,maxm,dx,dy,q_hbox_u1,aux_hbox_u1,alpha)
    call up_hboxes(x_0,y_0,i_0,j_0,x_e,y_e,i_e,j_e,2,&
        qold,aux,num_eqn,1,num_ghost,mx,my,maxm,dx,dy,q_hbox_u2,aux_hbox_u2,alpha)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (index_capa == 0) then
    !        # no capa array:
        do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
        5 END DO
    endif


    if( ids == 1 )then

    !     # perform x-sweeps
    !     ==================

    !     # note that for dimensional splitting we sweep over the rows of
    !     # ghosts cells as well as the interior.  This updates the ghost
    !     # cell values to the intermediate state as needed in the following
    !     # sweep in the y-direction.

        do 50 j = 1-num_ghost,my+num_ghost

        !        # copy data along a slice into 1d arrays:
            forall (m=1:num_eqn, i = 1-num_ghost: mx+num_ghost)
            q1d(m,i) = qold(m,i,j)
            end forall

            if (index_capa > 0)  then
                do 22 i = 1-num_ghost, mx+num_ghost
                    dtdx1d(i) = dtdx / aux(index_capa,i,j)
                22 END DO
            endif

            if (num_aux > 0)  then
                do 23 ma=1,num_aux
                    do 23 i = 1-num_ghost, mx+num_ghost
                        aux2(ma,i) = aux(ma,i,j  )
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
            q1d,dtdx1d,aux1,aux2,aux3,method,mthlim, &
            qadd,fadd,gadd,cfl1d, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave,vert,horz,bar_index_i,bar_ht,bar_loc,alpha)
            cfl = dmax1(cfl,cfl1d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if vertical barrier, then update hboxes here
            amdq_hbox = 0.d0
            apdq_hbox = 0.d0

            if (vert) then
              amdq_hbox(:,2) = amdq(:,bar_index+1) + f(q1d(:,bar_index),1)-f(q_hbox_d1(:,j),1)
              apdq_hbox(:,3) = apdq(:,bar_index+1) + f(q_hbox_u1(:,j),1)-f(q1d(:,bar_index+1),1)

              call riemann_solver(q_hbox_u1(:,j),q_hbox_u2(:,j),aux_hbox_u1(1,j),aux_hbox_u2(1,j),&
              s_r,fwave_r,amdq_hbox(:,3),apdq_hbox(:,4),1)

              amdq_hbox(:,4) = alpha*(f(q1d(:,bar_index+3),1)+amdq(:,bar_index+4)) &
              + (1-alpha)*(f(q1d(:,bar_index+2),1)+amdq(:,bar_index+3)) - f(q_hbox_u2(:,j),1)

              apdq_hbox(:,1) = f(q_hbox_d2(:,j),1) - &
              (alpha*(f(q1d(:,bar_index-1),1)-apdq(:,bar_index-1)) + (1-alpha)*(f(q1d(:,bar_index-2),1)-apdq(:,bar_index-2)))

              call riemann_solver(q_hbox_d2(:,j),q_hbox_d1(:,j),aux_hbox_d2(1,j),aux_hbox_d1(1,j),&
              s_r,fwave_r,amdq_hbox(:,1),apdq_hbox(:,2),1)
              ! update hboxes:
              q_hbox_d2(:,j)=q_hbox_d2(:,j)-dtdx*(apdq_hbox(:,1)+amdq_hbox(:,1))
              q_hbox_d1(:,j)=q_hbox_d1(:,j)-dtdx*(apdq_hbox(:,2)+amdq_hbox(:,2))
              q_hbox_u1(:,j)=q_hbox_u1(:,j)-dtdx*(apdq_hbox(:,3)+amdq_hbox(:,3))
              q_hbox_u2(:,j)=q_hbox_u2(:,j)-dtdx*(apdq_hbox(:,4)+amdq_hbox(:,4))
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !        # update qnew by flux differencing.
        !        # (rather than maintaining arrays f and g for the total fluxes,
        !        # the modifications are used immediately to update qnew
        !        # in order to save storage.)

            if (index_capa == 0) then

            !            # no capa array.  Standard flux differencing:
                forall (m=1:num_eqn, i=1:mx)
                    qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                    - dtdx * (fadd(m,i+1) - fadd(m,i))

                end forall
            else

            !            # with capa array.
                forall (m=1:num_eqn, i=1:mx)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                / aux(index_capa,i,j)
                end forall
            endif
        50 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! update q according to hboxes here:
          if (vert) then
          qnew(:,bar_index,:) = q_hbox_d1
          qnew(:,bar_index+1,:) = q_hbox_u1
          qnew(:,bar_index-1,:) = (1-alpha)*q_hbox_d1 + alpha*(q_hbox_d2)
          qnew(:,bar_index+2,:) = alpha*(q_hbox_u1) + (1-alpha)*q_hbox_u2
          qnew(:,bar_index-2,:) = temp
          qnew(:,bar_index+3,:) = temp
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif

    if( ids == 2 )then

    !     # perform y sweeps
    !     ==================


        do 100 i = 1-num_ghost, mx+num_ghost

        !        # copy data along a slice into 1d arrays:
            forall (m=1:num_eqn, j = 1-num_ghost: my+num_ghost)
            q1d(m,j) = qold(m,i,j)
            end forall

            if (index_capa > 0)  then
                do 72 j = 1-num_ghost, my+num_ghost
                    dtdy1d(j) = dtdy / aux(index_capa,i,j)
                72 END DO
            endif

            if (num_aux > 0)  then

                do 73 ma=1,num_aux
                    do 73 j = 1-num_ghost, my+num_ghost
                        aux2(ma,j) = aux(ma,i,j)
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
            q1d,dtdy1d,aux1,aux2,aux3,method,mthlim, &
            qadd,fadd,gadd,cfl1d, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave,vert,horz,bar_index_i,bar_ht,bar_loc,alpha)

            cfl = dmax1(cfl,cfl1d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if vertical barrier, then update hboxes here
            amdq_hbox = 0.d0
            apdq_hbox = 0.d0

            if (horz) then
              ! herin lies the seg fault!!!!!! the i+2 guys
              amdq_hbox(:,2) = amdq(:,bar_index+1) + f(q1d(:,bar_index),2)-f(q_hbox_d1(:,i),2)
              apdq_hbox(:,3) = apdq(:,bar_index+1) + f(q_hbox_u1(:,i),2)-f(q1d(:,bar_index+1),2)

              call riemann_solver(q_hbox_u1(:,i),q_hbox_u2(:,i),aux_hbox_u1(1,i),aux_hbox_u2(1,i),&
              s_r,fwave_r,amdq_hbox(:,3),apdq_hbox(:,4),2)
              amdq_hbox(:,4) = alpha*(f(q1d(:,bar_index+3),2)+amdq(:,bar_index+4)) &
              + (1-alpha)*(f(q1d(:,bar_index+2),2)+amdq(:,bar_index+3)) - f(q_hbox_u2(:,i),2)

              apdq_hbox(:,1) = f(q_hbox_d2(:,i),2) - &
              (alpha*(f(q1d(:,bar_index-1),2)-apdq(:,bar_index-1)) + (1-alpha)*(f(q1d(:,bar_index-2),2)-apdq(:,bar_index-2)))
              call riemann_solver(q_hbox_d2(:,i),q_hbox_d1(:,i),aux_hbox_d2(1,i),aux_hbox_d1(1,i),&
              s_r,fwave_r,amdq_hbox(:,1),apdq_hbox(:,2),2)
              ! update hboxes:
              q_hbox_d2(:,i)=q_hbox_d2(:,i)-dtdy*(apdq_hbox(:,1)+amdq_hbox(:,1))
              q_hbox_d1(:,i)=q_hbox_d1(:,i)-dtdy*(apdq_hbox(:,2)+amdq_hbox(:,2))
              q_hbox_u1(:,i)=q_hbox_u1(:,i)-dtdy*(apdq_hbox(:,3)+amdq_hbox(:,3))
              q_hbox_u2(:,i)=q_hbox_u2(:,i)-dtdy*(apdq_hbox(:,4)+amdq_hbox(:,4))
              ! write(*,*) "HBOX fluc: ", apdq_hbox(:,1),amdq_hbox(:,1)
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !        # update qnew by flux differencing.
        !        # Note that the roles of fadd and gadd are reversed for
        !        # the y-sweeps -- fadd is the modification to g-fluxes and
        !        # gadd is the modification to f-fluxes to the left and right.

            if (index_capa == 0) then
              forall (m=1:num_eqn,j=1:my)

                  qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                  - dtdx * (fadd(m,i+1) - fadd(m,i))
                end forall
            else

            !            # with capa array.
                forall (m=1:num_eqn, j=1:my)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                / aux(index_capa,i,j)
                end forall

            endif


        100 END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  update q according to hboxes here:
          if (horz) then
          qnew(:,:,bar_index) = q_hbox_d1
          qnew(:,:,bar_index+1) = q_hbox_u1
          qnew(:,:,bar_index-1) = (1-alpha)*q_hbox_d1 + alpha*(q_hbox_d2)
          qnew(:,:,bar_index+2) = alpha*(q_hbox_u1) + (1-alpha)*q_hbox_u2
          qnew(:,:,bar_index-2) = alpha*(qnew(:,:,bar_index-2)) + (1-alpha)*q_hbox_d2
          qnew(:,:,bar_index+3) = alpha*(q_hbox_u2) + (1-alpha)*qnew(:,:,bar_index+3)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif

    return
  end subroutine step2ds

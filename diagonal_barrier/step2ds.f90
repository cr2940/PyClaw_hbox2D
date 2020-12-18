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

    dimension amdqd2(num_eqn, mx)
    dimension waved2(num_eqn, num_waves, mx)
    dimension sd2(num_waves,mx)
    dimension apdqd2(num_eqn, mx)
    dimension amdqd1(num_eqn, mx+2)
    dimension waved1(num_eqn, num_waves, mx)
    dimension sd1(num_waves,mx)
    dimension apdqd1(num_eqn, mx+2)
    dimension amdqu1(num_eqn, mx+2)
    dimension waveu1(num_eqn, num_waves, mx)
    dimension su1(num_waves,mx)
    dimension apdqu1(num_eqn, mx+2)
    dimension amdqu2(num_eqn, mx)
    dimension waveu2(num_eqn, num_waves, mx)
    dimension su2(num_waves,mx)
    dimension apdqu2(num_eqn, mx)
    dimension s(3),fwave(3,3)
    dimension fwave_t(num_eqn,num_waves,mx+2), s_t(num_waves,mx+2)
    double precision :: amdq_hbox(num_eqn),apdq_hbox(num_eqn)
    double precision :: amdq_hbox2(num_eqn),apdq_hbox2(num_eqn)
    double precision :: tres_d1(3), tres_d2(3), tres_u1(3),tres_u2(3)
    double precision :: nout_d(3,mx+2), nout_u(3,mx+2)
    real(kind=8):: q_hbox_d1(3,mx+2),q_hbox_d2(3,mx+2)
    real(8) :: q_hbox_u2(3,mx+2)
    real(kind=8):: aux_hbox_d1(mx+2),q_hbox_u1(3,mx+2)
    real(8)::aux_hbox_u1(mx+2)
    real(kind=8):: aux_hbox_d2(mx+2), aux_hbox_u2(mx+2)
    real(8) :: H_12L(3),G_12L(3),H_12U(3),G_12U(3),E_L(3),E_U(3) ! outer fluxes of Lower hbox grid and Upper hbox grid
    real(8) :: H_32L(3),G_32L(3),H_32U(3),G_32U(3) ! outer fluxes (on the top part of grid, the lateral flux to last hboxes)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! barrier parameters
    real(kind=8) :: temp(3),temp2(3)
    real(kind=8)::x_0,x_e,y_0,y_e,alpha,bar_height,bar_ht,bar_loc
    integer :: i_0,i_e,j_0,j_e,bar_index,bar_index_i,k
    integer :: ind(mx+2),ind2(mx+2)
    logical :: L2R,R2L,xor_lr,check_on
    real(kind=8) :: xe(1-num_ghost:mx+num_ghost+1), ye(1-num_ghost:my+num_ghost+1)
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

    index_capa = method(6)
    num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! grid and barrier info
    bar_height = bar_ht
    ! print*,"bar ht:",bar_ht,bar_height
    ind = (/(i,i=1,mx+2)/)
    ind2 = ind + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get your hboxes & rotate  ! normalize by hbox area, sqrt(2)
    q_hbox_d1 = 0.d0
    q_hbox_u1 = 0.d0
    tres_d1 = 0.d0
    tres_u1 = 0.d0
    nout_d = 0.d0
    nout_u = 0.d0
    ! fragments weights:
      sqrt2 = sqrt(2.d0)
      a1 = 0.5d0  !alpha
      a2 = (1.d0-sqrt2/2.d0)*2 ! gamma
      a3= (sqrt2-(a1+a2))/(2.d0)! beta
      a4 = (0.5d0-a3) ! delta
      a5 = 1-a2-(2-1.5d0*sqrt2)*2
      a6 = (sqrt2 - a5 - 2.d0*a4)/(2.d0)! phi

      do i=1,mx+2  ! starts from second layer of ghost cell on left_bottom to first layer of ghost cell on right_top
      q_hbox_d1(:,i)=a1*qold(:,i-1,i-1)+a3*qold(:,i-1,i-2) &
     + a3*qold(:,i,i-1) + a2*qold(:,i,i-2)
       q_hbox_d1(:,i) = q_hbox_d1(:,i)/sqrt2
     q_hbox_u1(:,i)=a1*qold2(:,i-1,i-1)+a3*qold2(:,i-1,i)&
      + a3*qold2(:,i-2,i-1)+a2*qold2(:,i-2,i)
       q_hbox_u1(:,i) = q_hbox_u1(:,i)/sqrt2
      ! need to switch orientation if bouncing off the wall and not overtopping
      call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),&
       (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
      call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i), &
       (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))

    ! bathymetry averages of hboxes:
     aux_hbox_d1(i)=a1*aux(1,i-1,i-1)+a3*aux(1,i-1,i-2)&
      + a3*aux(1,i,i-1)+a2*aux(1,i,i-2)
      aux_hbox_d1(i) = aux_hbox_d1(i)/sqrt2
     aux_hbox_u1(i)=a1*auxu(1,i-1,i-1)+a3*auxu(1,i-1,i)&
      + a3*auxu(1,i-2,i-1)+a2*auxu(1,i-2,i)
      aux_hbox_u1(i) = aux_hbox_u1(i)/sqrt2
      end do

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !        # update qnew by flux differencing.
        !        # (rather than maintaining arrays f and g for the total fluxes,
        !        # the modifications are used immediately to update qnew
        !        # in order to save storage.)
        ! if (j .eq. 1-num_ghost) then
        !   tres_d1 =  a3*qadd(:,2-num_ghost)+a2*qadd(:,1)
        !   tres_d1 = 2.d0*tres_d1
        ! end if
        ! if (j .eq. my+num_ghost) then
        !   tres_u2 =  -a2*qadd2(:,mx) -a3*qadd2(:,mx+1)
        !   tres_u2 = 2.d0*tres_u2
        ! end if
        if ( j > 1-num_ghost .and. j < my+1) then
          if (j.ne.my) then
            nout_d(:,j+1) = a3*qadd(:,j+1)
            nout_d(:,j+2) = nout_d(:,j+2)+a2*qadd(:,j+2) + a3*qadd(:,j+1)
          else
            nout_d(:,j+1) =  a3*qadd(:,j+1)
          end if
        end if
        if ( j > 1-num_ghost .and. j < my+1) then
          if (j.ne.my) then
            nout_u(:,j+1) = a3*qadd2(:,j-1)
            nout_u(:,j+2) = nout_u(:,j+2)+a2*qadd2(:,j-2) + a3*qadd2(:,j-1)
          else
            nout_u(:,j+1) =  a3*qadd2(:,j-1)
          end if
        end if
        ! if ( j > 2-num_ghost) then
        !   if (j.ne.my+2) then
        !     nout_u(:,j) =  a2*qadd2(:,j-2)!+a3*qadd2(:,j-1)
        !     nout_u(:,j+1) = 2.d0*a3*qadd2(:,j-1)
        !   else
        !     nout_u(:,j) =  a2*qadd2(:,j-2)
        !   end if
        ! end if
        if (j .eq. 1-num_ghost) then
          tres_d1 =  a3*qadd(:,2-num_ghost)+a2*qadd(:,1)
          ! nout_d(:,1) = nout_d(:,1) + a2*qadd(:,1)
          tres_d1 = 2.d0*tres_d1
        end if
        if (j.eq.mx+1) then
          tres_d2 = a3*qadd(:,mx+num_ghost)
          tres_d2 = 2.d0*tres_d2
        end if
        if (j .eq. my+num_ghost) then
          tres_u2 = a3*qadd2(:,mx-1)
          tres_u2 = 2.d0*tres_u2
        end if
        if (j.eq.0) then
          tres_u1 =  a3*qadd2(:,-1)
          tres_u1 = 2.d0*tres_u1
        end if
        ! if ( j > 1-num_ghost .and. j < my+1) then
        !   if (j.eq.my) then
        !     nout_d(:,j+2) = nout_d(:,j+2) + a3*qadd(:,j+1)
        !     tres_d2 = tres_d2 + 2.d0*a2*qadd(:,j+2)
        !     nout_d(:,j+1) = nout_d(:,j+1) + a3*qadd(:,j+1)
        !   else
        !     nout_d(:,j+1) = nout_d(:,j+1) + a3*qadd(:,j+1)
        !     nout_d(:,j+2) = nout_d(:,j+2) + a3*qadd(:,j+1)+a2*qadd(:,j+2)
        !   end if
        ! end if
        ! if ( j > 2-num_ghost .and. j < mx+2) then
        !      if (j.eq.3-num_ghost) then
        !        tres_u1 = tres_u1 + 2.d0*a2*qadd2(:,j-2)
        !        nout_u(:,j) = nout_u(:,j) +a3*qadd2(:,j-1)
        !        nout_u(:,j+1) = nout_u(:,j+1) + a3*qadd2(:,j-1)
        !      else
        !       nout_u(:,j+1) = nout_u(:,j+1) +a3*qadd2(:,j-1)
        !       nout_u(:,j) = nout_u(:,j) +a2*qadd2(:,j-2)+a3*qadd2(:,j-1)
        !     end if
        ! end if


            if (index_capa == 0) then

            !            # no capa array.  Standard flux differencing:
                ! forall (m=1:num_eqn, i=1:mx)
                do m=1,num_eqn
                  do i=1,mx
                      qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                      - dtdx * (fadd(m,i+1) - fadd(m,i))
                      qnew2(m,i,j) = qnew2(m,i,j) + qadd2(m,i) &
                      - dtdx * (fadd2(m,i+1) - fadd2(m,i))
                  end do
                end do
                ! end forall
            else

            !            # with capa array.
                forall (m=1:num_eqn, i=1:mx)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                / auxu(index_capa,i,j)
                end forall
            endif
        50 END DO
    endif

    if( ids == 2 )then

    !     # perform y sweeps
    !     ==================

        do 100 i = 1-num_ghost, mx+num_ghost

        !        # copy data along a slice into 1d arrays:
            forall (m=1:num_eqn, j = 1-num_ghost: my+num_ghost)
            q1d(m,j) = qold(m,i,j)
            q1d2(m,j) = qold2(m,i,j)
            end forall
            if (index_capa > 0)  then
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
            call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my, &
            q1d2,dtdy1d,aux1,aux2u,aux3,method,mthlim, &
            qadd2,fadd2,gadd2,cfl1d2, &
            work(i0wave),work(i0s),amdq,apdq, &
            work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
            use_fwave)
            gm2(:,i,:) = amdq
            gp2(:,i,:) = apdq

            cfl = dmax1(cfl,min(cfl1d,cfl1d2))
        !        # update qnew by flux differencing.
        !        # Note that the roles of fadd and gadd are reversed for
        !        # the y-sweeps -- fadd is the modification to g-fluxes and
        !        # gadd is the modification to f-fluxes to the left and right.

          ! if (i.eq.mx+num_ghost) then
          !   tres_d2 = -a2*qadd(:,mx+num_ghost-2)-a3*qadd(:,mx+num_ghost-1)
          !   tres_d2 = 2.d0*tres_d2
          ! end if
          ! if (i.eq.-1) then
          !   tres_u1 =  a3*qadd2(:,0)+a2*qadd2(:,1)
          !   tres_u1 = 2.d0*tres_u1
          ! end if
          if ( i > 1-num_ghost .and. i < mx+1) then
            if (i.ne.mx) then
              nout_u(:,i+1) =  a3*qadd2(:,i+1)
              nout_u(:,i+2) =  nout_u(:,i+2)+a2*qadd2(:,i+2)+a3*qadd2(:,i+1)
            else
              nout_u(:,i+1) = nout_u(:,i+1)+a3*qadd2(:,i+1)
            end if
          end if
          if ( i > 1-num_ghost .and. i < mx+1) then
            if (i.ne.mx) then
              nout_d(:,i+1) = a3*qadd(:,i-1)
              nout_d(:,i+2) =  nout_d(:,i+2)+a2*qadd(:,i-2)+a3*qadd2(:,i-1)
            else
              nout_d(:,i+1) = nout_d(:,i+1)+a3*qadd(:,i-1)
            end if
          end if
          ! if ( i > 2-num_ghost ) then
          !   if (i.ne.mx+2) then
          !     nout_d(:,i) =  a2*qadd(:,i-2)!+a3*qadd(:,i-1)
          !     nout_d(:,i+1) = 2.d0*a3*qadd(:,i-1)
          !   else
          !     nout_d(:,i) =  a2*qadd(:,i-2)
          !   end if
          ! end if

          if (i.eq.mx+num_ghost) then
            tres_d2 = a3*qadd(:,mx+num_ghost-1)
            tres_d2 = 2.d0*tres_d2
          end if
          if (i.eq.-1) then
            tres_u1 =  a3*qadd2(:,0)+a2*qadd2(:,1)
            ! nout_u(:,1) = nout_u(:,1)+a2*qadd2(:,1)
            tres_u1 = 2.d0*tres_u1
          end if
          if (i .eq. 2-num_ghost) then
            tres_d1 = tres_d1 +  2*a3*qadd(:,1-num_ghost)
            ! tres_d1 = 2.d0*tres_d1
          end if
          if (i .eq. my+1) then
            tres_u2 =  a3*qadd2(:,my+2)
            tres_u2 = 2.d0*tres_u2
          end if
          ! if ( i > 1-num_ghost .and. i < mx+1) then
          !   if (i.eq.mx) then
          !     tres_u2 = tres_u2 + 2.d0*a2*qadd2(:,i+2)
          !     nout_u(:,i+2) = nout_u(:,i+2) + a3*qadd2(:,i+1)
          !     nout_u(:,i+1) = nout_u(:,i+1) + a3*qadd2(:,i+1)
          !   else
          !     nout_u(:,i+1) = nout_u(:,i+1) +a3*qadd2(:,i+1)
          !     nout_u(:,i+2) = nout_u(:,i+2) +a2*qadd2(:,i+2)+a3*qadd2(:,i+1)
          !   end if
          ! end if
          ! if ( i > 2-num_ghost .and. j < my+2) then
          !     if (i.eq.1) then
          !       nout_d(:,i) = nout_d(:,i) + a3*qadd(:,i-1)
          !       tres_d1 = tres_d1 + 2.d0*a2*qadd(:,i-2)
          !       nout_d(:,i+1) = nout_d(:,i+1) + a3*qadd(:,i-1)
          !     else
          !       nout_d(:,i+1) = nout_d(:,i+1) + a3*qadd(:,i-1)
          !       nout_d(:,i) = nout_d(:,i) + a3*qadd(:,i-1)+a2*qadd(:,i-2)
          !     end if
          ! end if


            if (index_capa == 0) then
              ! forall (m=1:num_eqn,j=1:my)
                do m=1,num_eqn
                  do j=1,my
                      qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                      - dtdx * (fadd(m,i+1) - fadd(m,i))
                      qnew2(m,i,j) = qnew2(m,i,j) + qadd2(m,j) &
                      - dtdx * (fadd2(m,i+1) - fadd2(m,i))
                  end do
                end do
                ! end forall

            else

            !            # with capa array.
                forall (m=1:num_eqn, j=1:my)
                qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                / auxu(index_capa,i,j)
                end forall

            endif


        100 END DO

    endif
    ! print* ,"FLUCTUATIONS: fp1>>", fp1
    ! print* ,"FLUCTUATIONS: fm1>>", fm1
    ! print* ,"FLUCTUATIONS: gp1>>", gp1
    ! print* ,"FLUCTUATIONS: gm1>>", gm1
    ! print* ,"FLUCTUATIONS: fp2>>", fp2
    ! print* ,"FLUCTUATIONS: fm2>>", fm2
    ! print* ,"FLUCTUATIONS: gp2>>", gp2
    ! print* ,"FLUCTUATIONS: gm2>>", gm2
      nout_u = 2.d0*nout_u
      nout_d = 2.d0*nout_d


    ! do i=1,mx+2
    !   call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),&
    !    (/-1/sqrt2,1/sqrt2/),(/-1/sqrt2,-1/sqrt2/))
    !   call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i), &
    !    (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
      ! call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i),(/-1/sqrt2,-1/sqrt2/),(/1/sqrt2,-1/sqrt2/))
      ! call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),(/1/sqrt2,1/sqrt2/),(/-1/sqrt2,1/sqrt2/))
    ! end do

    ! call rotate_state(tres_d1,tres_d1, &
    ! (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
    ! call rotate_state(tres_u1,tres_u1,&
    !  (/-1/sqrt2,1/sqrt2/),(/-1/sqrt2,-1/sqrt2/))
    !  call rotate_state(tres_d2,tres_d2, &
    !  (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
    !  call rotate_state(tres_u2,tres_u2,&
    !   (/-1/sqrt2,1/sqrt2/),(/-1/sqrt2,-1/sqrt2/))


     check_on = .false.
   ! normal:
    do i=1,mx+2
      hL = q_hbox_u1(1,i)
      hR = q_hbox_d1(1,i)
      huL = q_hbox_u1(2,i)
      huR = q_hbox_d1(2,i)
      hvL= q_hbox_u1(3,i)
      hvR = q_hbox_d1(3,i)
      bL = aux_hbox_u1(i)
      bR = aux_hbox_d1(i)
      temp = (/hL,huL,hvL/)!qold2(:,i-2,i-2)!
      temp2 = (/hR,huR,hvR/)!qold(:,i-2,i-2)!
      ! need to switch orientation if bouncing off the wall and not overtopping
      ! call rotate_state(temp,temp, &
      ! (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
      ! call rotate_state(temp2,temp2,&
      !  (/-1/sqrt2,1/sqrt2/),(/-1/sqrt2,-1/sqrt2/))
      call barrier_passing(temp(1),temp2(1),temp(2),temp2(2),&
                  bL,bR,bar_height,&
               L2R,R2L,hstarL,hstarR,ustarL,ustarR)
      xor_lr = XOR(L2R,R2L)
      amdq_hbox=0.d0
      apdq_hbox=0.d0
      amdq_hbox2=0.d0
      apdq_hbox2=0.d0

      ! need to switch orientation if bouncing off the wall and not overtopping
      ! if (xor_lr.or. (.not.L2R.and..not.R2L)) then
      ! call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i), &
      ! (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
      ! call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),&
      !  (/-1/sqrt2,1/sqrt2/),(/-1/sqrt2,-1/sqrt2/))
      !  end if
       ! check overtop in general
      if (L2R .or. R2L) then
        check_on = .true.
      end if
    call redistribute_fwave(1,temp,temp2,&
    aux_hbox_u1(i),aux_hbox_d1(i),bar_height,1,fwave,&
     s,amdq_hbox(:),apdq_hbox(:),num_waves,3,L2R,R2L)
     ! call redistribute_fwave(2,q_hbox_u1(:,i),q_hbox_d1(:,i),&
     ! aux_hbox_u1(i),aux_hbox_d1(i),bar_height,1,fwave,&
     !  s,amdq_hbox2(:),apdq_hbox2(:),num_waves,3,L2R,R2L)
     ! amdq_hbox = sqrt2*(amdq_hbox)!+amdq_hbox2)
     ! apdq_hbox = sqrt2*(apdq_hbox)!+apdq_hbox2)
     q_hbox_d1(:,i) = q_hbox_d1(:,i) - dtdx*apdq_hbox
     q_hbox_u1(:,i) = q_hbox_u1(:,i) - dtdx*amdq_hbox
     ! need to switch orientation if bouncing off the wall and not overtopping
     ! if (xor_lr.or. (.not.L2R.and..not.R2L)) then
     call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i),(/1/sqrt2,1/sqrt2/),(/-1/sqrt2,1/sqrt2/))
     call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),(/1/sqrt2,1/sqrt2/),(/-1/sqrt2,1/sqrt2/))
     ! end if
     q_hbox_d1(:,i) = q_hbox_d1(:,i) - nout_d(:,i)
     q_hbox_u1(:,i) = q_hbox_u1(:,i) + nout_u(:,i)
     end do


     do i=1,mx+2
       call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i), &
       (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
       call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i),&
        (/1/sqrt2,-1/sqrt2/),(/1/sqrt2,1/sqrt2/))
      end do
     ! transverse:
     call rn2(2, mx+2, num_eqn,num_waves,num_aux,num_ghost, &
              mx, q_hbox_d1, q_hbox_d1, aux_hbox_d1, aux_hbox_d1,&
               fwave_t, s_t, amdqd1, apdqd1)
     call rn2(2, mx+2, num_eqn,num_waves,num_aux,num_ghost, &
              mx, q_hbox_u1, q_hbox_u1, aux_hbox_u1, aux_hbox_u1,&
               fwave_t, s_t, amdqu1, apdqu1)
     forall(m=1:num_eqn,i=2:mx+2)
       q_hbox_d1(m,i) = q_hbox_d1(m,i) - dtdy*(1/sqrt2)*apdqd1(m,i)
       q_hbox_d1(m,i-1) = q_hbox_d1(m,i-1) - dtdy*(1/sqrt2)*amdqd1(m,i)
     end forall
     forall(m=1:num_eqn,i=2:mx+2)
       q_hbox_u1(m,i) = q_hbox_u1(m,i) - dtdy*(1/sqrt2)*apdqu1(m,i)
       q_hbox_u1(m,i-1) = q_hbox_u1(m,i-1) - dtdy*(1/sqrt2)*amdqu1(m,i)
     end forall
     do i=1,mx+2
       call rotate_state(q_hbox_u1(:,i),q_hbox_u1(:,i),(/1/sqrt2,1/sqrt2/),(/-1/sqrt2,1/sqrt2/))
       call rotate_state(q_hbox_d1(:,i),q_hbox_d1(:,i),(/1/sqrt2,1/sqrt2/),(/-1/sqrt2,1/sqrt2/))
     end do
     q_hbox_d1(:,1) = q_hbox_d1(:,1) + tres_d1
     q_hbox_d1(:,mx+2) = q_hbox_d1(:,mx+2) - tres_d2
     q_hbox_u1(:,1) = q_hbox_u1(:,1) + tres_u1
     q_hbox_u1(:,mx+2) = q_hbox_u1(:,mx+2) - tres_u2

     do i=1,mx+2
       qnew(:,i-1,i-1) = q_hbox_d1(:,i)!
       qnew2(:,i-1,i-1) = q_hbox_u1(:,i)
       qnew(:,i,i-2) = (1-a2)*qnew(:,i,i-2) + a2*q_hbox_d1(:,i)
       qnew2(:,i-2,i) = (1-a2)*qnew2(:,i-2,i) + a2*q_hbox_u1(:,i)
       if (i<mx+2.and.i>1) then
         qnew(:,i,i-1) = a3*q_hbox_d1(:,i) +a3*q_hbox_d1(:,i+1)+(1-2*a3)*qnew(:,i,i-1)
         qnew2(:,i-1,i) = a3*q_hbox_u1(:,i) +a3*q_hbox_u1(:,i+1)+(1-2*a3)*qnew2(:,i-1,i)
       end if
       if (i.eq.mx+2) then
         qnew(:,i,i-1) = a3*q_hbox_d1(:,i)+ (1-a3)*qnew(:,i,i-1)
         qnew2(:,i-1,i) = a3*q_hbox_u1(:,i) +(1-a3)*qnew2(:,i-1,i)
       end if
       if (i.eq.1) then
         qnew(:,i-1,i-2) = a3*q_hbox_d1(:,i)+(1-a3)*qnew(:,i-1,i-2)
         qnew2(:,i-2,i-1) =a3*q_hbox_u1(:,i)+(1-a3)*qnew2(:,i-2,i-1)
       end if
     end do
     ! qnew(:,-1,-1) = qnew(:,0,0)
     ! qnew2(:,-1,-1) = qnew2(:,0,0)
     ! qnew(:,mx+2,mx+2) = qnew(:,mx+1,mx+1)
     ! qnew2(:,mx+2,mx+2) = qnew2(:,mx+1,mx+1)

     ! ! ! exchange upper half values between qnew2 and qnew and lower half values between qnew and qnew2
    !  if (check_on) then
      ! do i=1,mx
      !   do j=1,my
      !     if (j<i) then
      !       qnew(:,i,j) = qnew2(:,i,j)
      !     else if (j>i) then
      !       qnew2(:,i,j) = qnew(:,i,j)
      !     end if
      !   end do
      ! end do
    ! else !if (.not.L2R .and. .not. R2L) then
      do i=1-num_ghost,mx+num_ghost
        do j=1-num_ghost,my+num_ghost
          if (j>i) then
            qnew(:,i,j) = qnew2(:,i,j)
          else if (j<i) then
            qnew2(:,i,j) = qnew(:,i,j)
          end if
        end do
      end do
    !  end if
     ! ! end if
       ! do i=1,mx
       !   do j=1,my
       !
       !     if (j>i) then
       !       print*,"Q_final:",q_final(1,:,:)
       !     end if
       !
       !   end do
       ! enddo
     !   do i=1,mx
     !     do j=1,my
     !       if (j<i) then
     !         qnew(:,i,j) = qnew2(:,i,j)
     !       else if (j>i) then
     !         qnew2(:,i,j) = qnew(:,i,j)
     !       end if
     !     end do
     !   end do
     ! end if

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

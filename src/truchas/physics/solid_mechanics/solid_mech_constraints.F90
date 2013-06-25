Module SOLID_MECH_CONSTRAINTS
  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  ! Contains:  
  !             
  ! Authors:  Dave Korzekwa (dak@lanl.gov)
  !-----------------------------------------------------------------------------
  Use kind_module,      Only: real_kind
  Use truchas_logging_services
!  Use parameter_module, Only: ndim
  Use var_vector_module
  use solid_mechanics_data
!  use nonlinear_solution, only: NK_SOLUTION_FIELD
  Implicit None

  Private

  ! Public procedures
  Public :: RHS_DISPLACEMENT_CONSTRAINTS,  &
            DISPLACEMENT_CONSTRAINTS,      &
            MECH_PRECOND_DISP_CONSTRAINTS, &
            FACE_GAP_INITIALIZE,           &
            FACE_GAP_UPDATE

  ! Penalty parameter for displacement BCs
  Real   (KIND = real_kind)                     :: penalty
  !-----------------------------------------------------------------------------

Contains

  Subroutine RHS_DISPLACEMENT_CONSTRAINTS()
    !=============================================================================
    !
    ! Apply nodal constraints such as interface and contact constraints,
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================

    use mech_bc_data_module
    use parameter_module,     only: ndim
    use bc_data_types

    implicit none
    ! Local variables
    Integer(KIND =  int_kind)                          :: inode, yindex, nnum, idim
    real(real_kind),dimension(ndim)                    :: Nvec1, Tvec, Avec, rj
    real(real_kind)                                    :: ndotr, tdotr, tdota, ndotdn, d

    penalty = contact_penalty

    ! Apply projections of source and enforced displacement vectors, making use of the 
    ! fact that the ith component of the vector [n n_T] u == n_i (n dot u)
    do inode = 1,SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       yindex = (nnum -1)*ndim
       select case(Node_Displacement_BC%Combination(inode))
       case(ONE_DISPLACEMENT, ONE_D_ONE_NC, ONE_D_TWO_NC)
          d = Node_Displacement_BC%Value(inode,1)
          ! c * d * [nn^T], expanded
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          rj = Src(yindex+1:yindex + ndim)
          ! Source terms
          ndotr = 0.0
          ! Displacement terms
          ndotdn = 0.0
          ! Accumulate n dot r and n dot dn
          ndotr = SUM(Nvec1(:) * rj(:))
          ndotdn = SUM(Nvec1(:) * d * Nvec1(:))
          ! Subtract source component normal to the surface
          do idim = 1,ndim
             RHS(yindex + idim) = -Src(yindex + idim) + ndotr * Nvec1(idim)
          end do
          ! Displacement terms
          do idim = 1,ndim
             RHS(yindex + idim) = RHS(yindex + idim) - penalty * ndotdn * Nvec1(idim)
          end do

       case(TWO_DISPLACEMENTS, TWO_D_ONE_NC)
          ! Combined displacement vector
          Avec = Node_Displacement_BC%Vector(inode,1,:)
          ! Tangent vector t \propto n1 x n2
          Tvec = Node_Displacement_BC%Vector(inode,2,:)
          ! Accumulate t dot r and t dot a
          rj = Src(yindex+1:yindex + ndim)
          tdotr = 0.0
          tdota = 0.0
          tdotr = SUM(Tvec(:) * rj(:))
          tdota = SUM(Tvec(:) * Avec(:))
          ! Start with tangent component of source terms
          do idim = 1,ndim
             RHS(yindex + idim) = - tdotr * Tvec(idim)
          end do
          ! Add in-plane component of a
          do idim = 1,ndim
             RHS(yindex + idim) = RHS(yindex + idim) - penalty * (Avec(idim) - tdota * Tvec(idim))
          end do

       case(THREE_DISPLACEMENTS)
          ! We already have a unique displacement vector in global coordinates
          do idim = 1,ndim
             RHS(yindex + idim) = penalty * Node_Displacement_BC%Vector(inode,1,idim)
          end do

       case(ONE_NORM_CONST)
          ! Nothing to do here: RHS = - Src

       case(TWO_NORM_CONST)
          ! Nothing to do here: RHS = - Src

       case(THREE_NORM_CONST)
          ! Nothing to do here: RHS = - Src

       end select
    end do

  end Subroutine RHS_DISPLACEMENT_CONSTRAINTS


  SUBROUTINE DISPLACEMENT_CONSTRAINTS(X,Y)
    !=============================================================================
    !
    ! Apply nodal constraints such as interface and contact constraints.
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================

    use mech_bc_data_module
    use parameter_module,     only: ndim, nnodes
    use bc_data_types
    use gs_module,            only: NN_Gather_BoundaryData
    use mesh_module,          only: Vertex

    implicit none

    ! Arguments
    Real   (KIND = real_kind), Dimension(:), Intent(IN)    :: X
    Real   (KIND = real_kind), Dimension(:), Intent(INOUT) :: Y

    ! Local variables
    integer                                    :: status
    Integer(KIND =  int_kind)                  :: inode, yindex, nnum, idim, gap_nnum1, &
                                                  gap_nnum2, gap_nnum3, gapindex1, gapindex2
    real(real_kind),dimension(ndim)            :: Nvec1, Tvec1, Tvec2, Tvec3, &
                                                  Gvec1, Gvec2, Gvec3, Vvec, Wvec, &
                                                  uj, uk, ul, fj, frk, frl, cj, ck, cl
    real(real_kind)                            :: ndotfj, ndotuj, ndotfrk, ndotukuj, ndotfrl, ndotuluj, &
                                                  tdotfj, tdotuj, tdotfrk, tdotfrl, tdotukuj, tdotuluj, &
                                                  t1dotfrk, t1dotukuj, &
                                                  t2dotfrk, t2dotukuj, &
                                                  t3dotfrk, t3dotukuj, & 
                                                  mdotfrk, mdotfrl, mdotukuj, mdotuluj, &
                                                  pdotfrk, pdotukuj, &
                                                  vdotfrk, vdotukuj, &
                                                  wdotfrk, wdotfrl, wdotukuj, wdotuluj, &
                                                  costheta1, costheta2, &
                                                  lambda1, lambda2, lambda3
    real(real_kind), dimension(ndim*nnodes)    :: Ysave, Ftot
    real(real_kind), pointer, dimension(:)     :: Btemp, X_Bound, XC_Bound, Y_Bound, S_Bound, CS_Bound
     
    NULLIFY(Btemp)
    NULLIFY(X_Bound, X_Bound, XC_Bound, Y_Bound, S_Bound, CS_Bound)

    penalty = contact_penalty

    ! Unmodified copy of Y
    Ysave = Y

    ! Gather displacement components, coordinates, matvec values and source terms (RHS) of 
    ! off-processor neighboring nodes
    do idim = 1,ndim
       call NN_GATHER_BOUNDARYDATA (SOURCE=X(idim:ndim*(nnodes-1)+idim:ndim),BOUNDARY=Btemp)
       if (idim == 1) then
          allocate(X_Bound(ndim * SIZE(Btemp)),stat = status)
          if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINTS: allocation error: X_Bound')
       end if
       X_Bound(idim:ndim*(SIZE(Btemp)-1)+idim:ndim) = Btemp
       DEALLOCATE(Btemp)
    end do

    do idim = 1,ndim
       call NN_GATHER_BOUNDARYDATA (SOURCE=Vertex%Coord(idim),BOUNDARY=Btemp)
       if (idim == 1) then
          allocate(XC_Bound(ndim * SIZE(Btemp)),stat = status)
          if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINTS: allocation error: XC_Bound')
       end if
       XC_Bound(idim:ndim*(SIZE(Btemp)-1)+idim:ndim) = Btemp
       DEALLOCATE(Btemp)
    end do


    do idim = 1,ndim
       call NN_GATHER_BOUNDARYDATA (SOURCE=Y(idim:ndim*(nnodes-1)+idim:ndim),BOUNDARY=Btemp)
       if (idim == 1) then
          allocate(Y_Bound(ndim * SIZE(Btemp)),stat = status)
          if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINTS: allocation error: Y_Bound')
       end if
       Y_Bound(idim:ndim*(SIZE(Btemp)-1)+idim:ndim) = Btemp
       DEALLOCATE(Btemp)
    end do

    do idim = 1,ndim
       call NN_GATHER_BOUNDARYDATA (SOURCE=Src(idim:ndim*(nnodes-1)+idim:ndim),BOUNDARY=Btemp)
       if (idim == 1) then
          allocate(S_Bound(ndim * SIZE(Btemp)),stat = status)
          if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINTS: allocation error: S_Bound')
       end if
       S_Bound(idim:ndim*(SIZE(Btemp)-1)+idim:ndim) = Btemp
       DEALLOCATE(Btemp)
    end do

    call NN_GATHER_BOUNDARYDATA (SOURCE=cscale(ndim:ndim*nnodes:ndim),BOUNDARY=Btemp)
    allocate(CS_Bound(SIZE(Btemp)),stat = status)
    if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINTS: allocation error: CS_Bound')
    CS_Bound = Btemp
    DEALLOCATE(Btemp)

    ! Get contact function values, lambda, for each gap node pair
    Ftot = Ysave + Src

    CALL CONTACT_FUNCTION(X, Ftot, X_Bound, XC_Bound)

    ! Apply projections of force and displacement vectors, making use of the 
    ! fact that the ith component of the vector [n n_T] u == n_i (n dot u)
    do inode = 1,SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       yindex = (nnum -1)*ndim
       select case(Node_Displacement_BC%Combination(inode))
       case(ONE_DISPLACEMENT)
          ! Normal vector
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Accumulate n dot f_j and n dot u_j
          ndotfj = 0.0
          ndotuj = 0.0
          uj = X(yindex+1:yindex+ndim)
          fj = Ysave(yindex+1:yindex + ndim)
          ndotfj = SUM(Nvec1(:) * fj(:))
          ndotuj = SUM(Nvec1(:) * uj(:))
          ! Subtract source component normal to the surface
          do idim = 1,ndim
             Y(yindex + idim) = Ysave(yindex + idim) - ndotfj * Nvec1(idim)
          end do
          ! Displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) - penalty * ndotuj * Nvec1(idim)
          end do

       case(TWO_DISPLACEMENTS)
          ! Tangent vector t \propto n1 x n2
          Tvec1 = Node_Displacement_BC%Vector(inode,2,:)
          ! Accumulate t dot f_j and t dot u_j
          tdotfj = 0.0
          tdotuj = 0.0
          uj = X(yindex+1:yindex+ndim)
          fj = Ysave(yindex+1:yindex + ndim)
          tdotfj = SUM(Tvec1(:) * fj(:))
          tdotuj = SUM(Tvec1(:) * uj(:))
          ! Start with tangent component of f_j
          do idim = 1,ndim
             Y(yindex + idim) = tdotfj * Tvec1(idim)
          end do
          ! Add in-plane component of u
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) - &
                  penalty * (X(yindex + idim) - tdotuj * Tvec1(idim))
          end do

       case(THREE_DISPLACEMENTS)
          ! We already have a unique displacement vector in global coordinates
          ! The displacement values are set in the right hand side
          do idim = 1,ndim
             Y(yindex + idim) = penalty * X(yindex + idim)
          end do

       case(ONE_NORM_CONST)
          ! Gap node number and gap index
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,1)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          ! Normal vector
          Gvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Need to initialize so that some force terms will be ignored
          fj = 0.0
          frk = 0.0
          uj = 0.0
          uk = 0.0
          cj = 0.0
          ck = 0.0
          ! Displacement and force vectors 
          uj = X(yindex+1:yindex + ndim)
          cj = Vertex(nnum)%Coord
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if

          lambda1 = Node_Displacement_BC%Lambda(inode,1)
          ! Accumulate n dot (f_k + r_k) and n dot (u_k - u_j)
          ndotfrk = 0.0
          ndotukuj = 0.0
          ndotfrk = SUM(Gvec1(:) * frk(:))
          ndotukuj = SUM(Gvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          ! Add force correction from Y and Source terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * (Gvec1(idim) * ndotfrk)
          end do
          ! Add displacement portion of force correction
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * penalty * (Gvec1(idim) * ndotukuj)
          end do

       case(TWO_NORM_CONST)
          ! Gap node numbers and gap indexes
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,1)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          gap_nnum2 = Node_Displacement_BC%Gap_Node(inode,2)
          if (gap_nnum1 /= gap_nnum2) then
             gapindex2 = (abs(gap_nnum2)-1)*ndim
          else
             ! Flag the condition where there is only one gap node
             gapindex2 = -1
          end if
          ! Gap normal vectors
          Gvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Gvec2 = Node_Displacement_BC%Normal(inode,2,:)
          ! Tangent vector
          Tvec1 = Node_Displacement_BC%Vector(inode,1,:)
          ! Need to initialize so that some force terms can be ignored
          fj = 0.0
          frk = 0.0
          frl = 0.0
          uj = 0.0
          uk = 0.0
          ul = 0.0
          cj = 0.0
          ck = 0.0
          cl = 0.0
          ! Displacement and force vectors 
          uj = X(yindex+1:yindex + ndim)
          cj = Vertex(nnum)%Coord
          fj = Ysave(yindex+1:yindex + ndim)
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if
          if (gapindex2 >= 0) then
             if (gap_nnum2 > 0) then
                ul = X(gapindex2+1:gapindex2+ndim)
                cl = Vertex(gap_nnum2)%Coord
                frl = (Ysave(gapindex2+1:gapindex2+ndim) + Src(gapindex2+1:gapindex2+ndim)) * &
                     cscale(gap_nnum2*ndim)/cscale(nnum * ndim)
             else
                ul = X_bound(gapindex2+1:gapindex2+ndim)
                cl = XC_bound(gapindex2+1:gapindex2+ndim)
                frl = (Y_Bound(gapindex2+1:gapindex2+ndim) + S_Bound(gapindex2+1:gapindex2+ndim)) * &
                     CS_Bound(abs(gap_nnum2))/cscale(nnum * ndim)
             end if
          end if
          !
          lambda1 = Node_Displacement_BC%Lambda(inode,1)
          lambda2 = Node_Displacement_BC%Lambda(inode,2)
          ! lambda3 = 1 unless Tvec = 0
          lambda3 = 1.0
!  All other use of lambda3 is bogus (?)
!          if (gapindex2 >= 0) lambda3 = Node_Displacement_BC%Lambda(inode,3)

          ! If the two normal vectors are the same, then Tvec = [0,0,0]. In this case 
          ! we set lambda3 to 0.
          if (ALL(Tvec1 == 0.0))  lambda3 = 0.0

          ! Initialize dot products to zero.  The logic for handling one or two gap nodes
          ! depends on this.
          ndotfrk = 0.0
          ndotukuj = 0.0
          mdotfrk = 0.0
          mdotfrl = 0.0
          mdotukuj = 0.0
          mdotuluj = 0.0
          tdotfrk = 0.0
          tdotfrl = 0.0
          tdotukuj = 0.0
          tdotuluj = 0.0
          ! Dot products
          ndotfrk = SUM(Gvec1(:) * frk(:))
          ndotukuj = SUM(Gvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          tdotfrk = SUM(Tvec1(:) * frk(:))
          tdotukuj = SUM(Tvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          if (gapindex2 < 0) then
             mdotfrk = SUM(Gvec2(:) * frk(:))
             mdotukuj = SUM(Gvec2(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          else
             mdotfrl = SUM(Gvec2(:) * frl(:))
             mdotuluj = SUM(Gvec2(:) * (ul(:) + cl(:) - uj(:) - cj(:)))
             tdotfrl = SUM(Tvec1(:) * frl(:))
             tdotuluj = SUM(Tvec1(:) * (ul(:) + cl(:) - uj(:) - cj(:)))
          end if

          ! Add force correction from Y and Source terms.  Note that mdotfrl and tdotfrl are zero 
          ! if there is only one gap node, and mdotfrk is zero if there are two gap nodes.
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * Gvec1(idim) * ndotfrk + &
                                lambda2 * Gvec2(idim) * (mdotfrk + mdotfrl) + &
                                lambda1 * lambda2 * lambda3 * (frk(idim) + frl(idim) - &
                                Tvec1(idim) * (tdotfrk + tdotfrl) - Gvec1(idim) * ndotfrk - &
                                Gvec2(idim) * (mdotfrk + mdotfrl))
          end do
          ! Add displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * penalty * Gvec1(idim) * ndotukuj + &
                                lambda2 * penalty * Gvec2(idim) * (mdotukuj + mdotuluj) + &
                                lambda1 * lambda2 * lambda3 * penalty * ((uk(idim) + ck(idim) - uj(idim) - cj(idim)) - &
                                Tvec1(idim) * (tdotukuj + tdotuluj) - Gvec1(idim) * ndotukuj - &
                                Gvec2(idim) * (mdotukuj + mdotuluj))
          end do
          ! Add displacement term for second gap node if it exists
          if (gapindex2 >= 0) then
             do idim = 1,ndim
                Y(yindex + idim) = Y(yindex + idim) + &
                                   lambda1 * lambda2 * lambda3 * penalty * (ul(idim) + cl(idim) - uj(idim) - cj(idim))
             end do
          end if


       case(THREE_NORM_CONST)
          ! Gap node numbers and gap indexes
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,1)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          gap_nnum2 = Node_Displacement_BC%Gap_Node(inode,2)
          gap_nnum3 = Node_Displacement_BC%Gap_Node(inode,3)
          if ((gap_nnum1 /= gap_nnum2) .or. (gap_nnum2 /= gap_nnum3)) then
             call TLS_panic ('DISPLACEMENT_CONSTRAINTS: We can only handle one gap node if there are three interfaces at a point')
          end if
          ! Gap normal vectors
          Gvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Gvec2 = Node_Displacement_BC%Normal(inode,2,:)
          Gvec3 = Node_Displacement_BC%Normal(inode,3,:)
          ! Tangent vectors
          Tvec1 = Node_Displacement_BC%Vector(inode,1,:)
          Tvec2 = Node_Displacement_BC%Vector(inode,2,:)
          Tvec3 = Node_Displacement_BC%Vector(inode,3,:)
          ! Need to initialize so that some force terms will be ignored
          fj = 0.0
          frk = 0.0
          uj = 0.0
          uk = 0.0
          ! Displacement and force vectors 
          uj = X(yindex+1:yindex + ndim)
          cj = Vertex(nnum)%Coord
          fj = Ysave(yindex+1:yindex + ndim)
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if
          lambda1 = Node_Displacement_BC%Lambda(inode,1)
          lambda2 = Node_Displacement_BC%Lambda(inode,2)
          lambda3 = Node_Displacement_BC%Lambda(inode,3)
          ndotfrk = 0.0
          ndotukuj = 0.0
          mdotfrk = 0.0
          mdotukuj = 0.0
          pdotfrk = 0.0
          pdotukuj = 0.0
          t1dotfrk = 0.0
          t1dotukuj = 0.0
          t2dotfrk = 0.0
          t2dotukuj = 0.0
          t3dotfrk = 0.0
          t3dotukuj = 0.0
          ! Dot products
          ndotfrk = SUM(Gvec1(:) * frk(:))
          ndotukuj = SUM(Gvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          mdotfrk = SUM(Gvec2(:) * frk(:))
          mdotukuj = SUM(Gvec2(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          pdotfrk = SUM(Gvec3(:) * frk(:))
          pdotukuj = SUM(Gvec3(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          t1dotfrk = SUM(Tvec1(:) * frk(:))
          t1dotukuj = SUM(Tvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          t2dotfrk = SUM(Tvec2(:) * frk(:))
          t2dotukuj = SUM(Tvec2(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          t3dotfrk = SUM(Tvec3(:) * frk(:))
          t3dotukuj = SUM(Tvec3(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          ! Force terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * Gvec1(idim) * ndotfrk + &
                                lambda2 * Gvec2(idim) * mdotfrk + &
                                lambda3 * Gvec3(idim) * pdotfrk + &
                                lambda1 * lambda2 * (frk(idim) - Tvec1(idim) * t1dotfrk - &
                                                 Gvec1(idim) * ndotfrk - Gvec2(idim) * mdotfrk) + &
                                lambda2 * lambda3 * (frk(idim) - Tvec2(idim) * t2dotfrk - &
                                                 Gvec2(idim) * mdotfrk - Gvec3(idim) * pdotfrk) + &
                                lambda3 * lambda1 * (frk(idim) - Tvec3(idim) * t3dotfrk - &
                                                 Gvec3(idim) * pdotfrk - Gvec1(idim) * ndotfrk) - &
                                lambda1 * lambda2 * lambda3 * (2.0 * frk(idim) - &
                                Tvec1(idim) * t1dotfrk - Tvec2(idim) * t2dotfrk - Tvec3(idim) * t3dotfrk - &
                                Gvec1(idim) * ndotfrk - Gvec2(idim) * mdotfrk - Gvec3(idim) * pdotfrk)
          end do
          ! Displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * penalty * Gvec1(idim) * ndotukuj + &
                                lambda2 * penalty * Gvec2(idim) * mdotukuj + &
                                lambda3 * penalty * Gvec3(idim) * pdotukuj + &
                                lambda1 * lambda2 * penalty * (uk(idim) + ck(idim) - uj(idim) - cj(idim) - Tvec1(idim) * &
                                    t1dotukuj -  Gvec1(idim) * ndotukuj - Gvec2(idim) * mdotukuj) + &
                                lambda2 * lambda3 * penalty * (uk(idim) + ck(idim) - uj(idim) - cj(idim) - Tvec2(idim) * &
                                    t2dotukuj - Gvec2(idim) * mdotukuj - Gvec3(idim) * pdotukuj) + &
                                lambda3 * lambda1 * penalty * (uk(idim) + ck(idim) - uj(idim) - cj(idim) - Tvec3(idim) * &
                                    t3dotukuj - Gvec3(idim) * pdotukuj - Gvec1(idim) * ndotukuj) - &
                                lambda1 * lambda2 * lambda3 * penalty * (2.0 * (uk(idim) + ck(idim) - uj(idim) - cj(idim)) - &
                                Tvec1(idim) * t1dotukuj - Tvec2(idim) * t2dotukuj - Tvec3(idim) * t3dotukuj - &
                                Gvec1(idim) * ndotukuj - Gvec2(idim) * mdotukuj - Gvec3(idim) * pdotukuj)
          end do

       case(ONE_D_ONE_NC)
          ! Gap node number and gap index
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,2)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          ! Normal vector
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Third ortho vector
          Vvec = Node_Displacement_BC%Vector(inode,2,:)
          ! Gap normal
          Gvec1 = Node_Displacement_BC%Normal(inode,2,:)
          ! Cosine theta (v dot n_g)
          costheta1 = Node_Displacement_BC%Scalar(inode,1)
          uj = X(yindex+1:yindex+ndim)
          cj = Vertex(nnum)%Coord
          fj = Ysave(yindex+1:yindex+ndim)
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_Bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if
          !
          lambda1 = Node_Displacement_BC%Lambda(inode,2)
          ! Accumulate n dot f_j, n dot u_j, v dot (f_k + r_k), v dot (u_k - u_j)
          ndotfj = 0.0
          ndotuj = 0.0
          vdotfrk = 0.0
          vdotukuj = 0.0

          ndotfj = SUM(Nvec1(:) * fj(:))
          ndotuj = SUM(Nvec1(:) * uj(:))
          vdotfrk = SUM(Vvec(:) * frk(:))
          vdotukuj = SUM(Vvec(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          ! Force terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) - Nvec1(idim) * ndotfj + lambda1 * vdotfrk * Vvec(idim)
          end do
          ! Displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + &
                                lambda1 * penalty * costheta1 * costheta1 * vdotukuj * Vvec(idim) - &
                                penalty * ndotuj * Nvec1(idim)
          end do

       case(TWO_D_ONE_NC)
          ! Gap node number and gap index
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,3)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          ! Tangent vector (displacement vector A is in the first index)
          Tvec1 = Node_Displacement_BC%Vector(inode,2,:)
          ! Gap normal
          Gvec1 = Node_Displacement_BC%Normal(inode,3,:)
          ! Cosine theta (t dot n_g)
          costheta1 = Node_Displacement_BC%Scalar(inode,1)
          uj = X(yindex+1:yindex + ndim)
          cj = Vertex(nnum)%Coord
          fj = Ysave(yindex+1:yindex + ndim)
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if
          !
          lambda1 = Node_Displacement_BC%Lambda(inode,3)
          ! Accumulate t dot f_j, t dot u_j, t dot (f_k + r_k) and t dot (u_k - u_j)
          tdotfj = 0.0
          tdotuj = 0.0
          tdotfrk = 0.0
          tdotukuj = 0.0
          tdotfj = SUM(Tvec1(:) * fj(:))
          tdotuj = SUM(Tvec1(:) * uj(:))
          tdotfrk = SUM(Tvec1(:) * frk(:))
          tdotukuj = SUM(Tvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          ! Force terms 
          do idim = 1,ndim
             Y(yindex + idim) = tdotfj * Tvec1(idim) + lambda1 * tdotfrk * Tvec1(idim)
          end do
          ! Displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + &
                                lambda1 * penalty * costheta1 * costheta1 * tdotukuj * Tvec1(idim) - &
                                penalty * (X(yindex + idim) - tdotuj * Tvec1(idim))
          end do

       case(ONE_D_TWO_NC) 
          ! Gap node numbers and gap indexes
          gap_nnum1 = Node_Displacement_BC%Gap_Node(inode,2)
          gapindex1 = (abs(gap_nnum1)-1)*ndim
          gap_nnum2 = Node_Displacement_BC%Gap_Node(inode,3)
          if (gap_nnum1 /= gap_nnum2) then
             gapindex2 = (abs(gap_nnum2)-1)*ndim
          else
             gapindex2 = -1
          end if
          ! Gap vectors in the plane of the free surface
          Vvec = Node_Displacement_BC%Vector(inode,1,:)
          Wvec = Node_Displacement_BC%Vector(inode,2,:)
          ! Surface normal vector
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Scalars
          costheta1 = Node_Displacement_BC%Scalar(inode,1)
          costheta2 = Node_Displacement_BC%Scalar(inode,2)
          ! Need to initialize so that some force terms will be ignored
          fj = 0.0
          frk = 0.0
          frl = 0.0
          uj = 0.0
          cj = 0.0
          uk = 0.0
          ck = 0.0
          ul = 0.0
          cl = 0.0
          ! Displacement and force vectors 
          uj = X(yindex+1:yindex + ndim)
          cj = Vertex(nnum)%Coord
          fj = Ysave(yindex+1:yindex + ndim)
          if (gap_nnum1 > 0) then
             uk = X(gapindex1+1:gapindex1+ndim)
             ck = Vertex(gap_nnum1)%Coord
             frk = (Ysave(gapindex1+1:gapindex1+ndim) + Src(gapindex1+1:gapindex1+ndim)) * &
                  cscale(gap_nnum1*ndim)/cscale(nnum * ndim)
          else
             uk = X_bound(gapindex1+1:gapindex1+ndim)
             ck = XC_bound(gapindex1+1:gapindex1+ndim)
             frk = (Y_Bound(gapindex1+1:gapindex1+ndim) + S_Bound(gapindex1+1:gapindex1+ndim)) * &
                  CS_Bound(abs(gap_nnum1))/cscale(nnum * ndim)
          end if
          if (gapindex2 >= 0) then
             if (gap_nnum2 > 0) then
                ul = X(gapindex2+1:gapindex2+ndim)
                cl = Vertex(gap_nnum2)%Coord
                frl = (Ysave(gapindex2+1:gapindex2+ndim) + Src(gapindex2+1:gapindex2+ndim)) * &
                     cscale(gap_nnum2*ndim)/cscale(nnum * ndim)
             else
                ul = X_bound(gapindex2+1:gapindex2+ndim)
                cl = XC_bound(gapindex2+1:gapindex1+ndim)
                frl = (Y_Bound(gapindex2+1:gapindex2+ndim) + S_Bound(gapindex2+1:gapindex2+ndim)) * &
                     CS_Bound(abs(gap_nnum2))/cscale(nnum * ndim)
             end if
          end if
          !
          lambda1 = Node_Displacement_BC%Lambda(inode,2)
          lambda2 = Node_Displacement_BC%Lambda(inode,3)
          ! If there is only one gap node lambda3 = 1
          lambda3 = 1.0
! lambda3 is always one unless the two normals are the same
          ! If the two normal vectors are the same, then Vvec = Wvec. In this case 
          ! we set lambda3 to 0.
          if (ALL(Vvec(:) == Wvec(:))) lambda3 = 0.0

          ! Initialize dot products to zero.  The logic for handling one or two gap nodes
          ! depends on this.
          ndotfj = 0.0
          ndotfrk = 0.0
          ndotfrl = 0.0
          ndotukuj = 0.0
          ndotuluj = 0.0
          vdotfrk = 0.0
          vdotukuj = 0.0
          wdotfrk = 0.0
          wdotfrl = 0.0
          wdotukuj = 0.0
          wdotuluj = 0.0
          ! Dot products
          ndotfj = SUM(Nvec1(:) * fj(:))
          ndotfrk = SUM(Nvec1(:) * frk(:))
          ndotuj = SUM(Nvec1(:) * uj(:))
          ndotukuj = SUM(Nvec1(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          vdotfrk = SUM(Vvec(:) * frk(:))
          vdotukuj = SUM(Vvec(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          if (gapindex2 < 0) then
             wdotfrk = SUM(Wvec(:) * frk(:))
             wdotukuj = SUM(Wvec(:) * (uk(:) + ck(:) - uj(:) - cj(:)))
          else
             ndotfrl = SUM(Nvec1(:) * frl(:))
             ndotuluj = SUM(Nvec1(:) * (ul(:) + cl(:) - uj(:) - cj(:)))
             wdotfrl = SUM(Wvec(:) * frl(:))
             wdotuluj = SUM(Wvec(:) * (ul(:) + cl(:) - uj(:) - cj(:)))
          end if
          ! Add force correction from Y and Source terms.  Note that ndotfrl and wdotfrl are zero 
          ! if there is only one gap node, and wdotfrk is zero if there are two gap nodes.
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) - Nvec1(idim) * ndotfj + lambda1 * Vvec(idim) * vdotfrk + &
                                lambda2 * Wvec(idim) * (wdotfrk + wdotfrl) + &
                                lambda1 * lambda2 * lambda3 * (frk(idim) + frl(idim) - &
                                Nvec1(idim) * (ndotfrk + ndotfrl) - Vvec(idim) * vdotfrk - &
                                Wvec(idim) * (wdotfrk + wdotfrl))
          end do
          ! Add displacement terms
          do idim = 1,ndim
             Y(yindex + idim) = Y(yindex + idim) + lambda1 * costheta1 * costheta1 * penalty * Vvec(idim) * vdotukuj + &
                                lambda2 * costheta2 * costheta2 * penalty * Wvec(idim) * (wdotukuj + wdotuluj) + &
                                lambda1 * lambda2 * lambda3 * penalty * (costheta1 * costheta1 * &
                                ((uk(idim) + ck(idim) - uj(idim) - cj(idim)) - &
                                Nvec1(idim) * (ndotukuj + ndotuluj) - Vvec(idim) * vdotukuj) - &
                                costheta2 * costheta2 * Wvec(idim) * (wdotukuj + wdotuluj)) - &
                                penalty * ndotuj * Nvec1(idim)
          end do
          ! Add term for second gap node if it exists
          if (gapindex2 >= 0) then
             do idim = 1,ndim
                Y(yindex + idim) = Y(yindex + idim) + &
                                   lambda1 * lambda2 * lambda3 * costheta2 * costheta2 * penalty * &
                                   (ul(idim) + cl(idim) - uj(idim) - cj(idim))
             end do
          end if

       end select
    end do

    DEALLOCATE(X_Bound)
    DEALLOCATE(XC_Bound)
    DEALLOCATE(Y_Bound)
    DEALLOCATE(S_Bound)
    DEALLOCATE(CS_Bound)

  END SUBROUTINE DISPLACEMENT_CONSTRAINTS

  SUBROUTINE CONTACT_FUNCTION(U, F, U_Bound, XC_Bound)

    use mech_bc_data_module
    use parameter_module,     only: ndim
    use mesh_module,          only: Vertex

    real(real_kind), dimension(:), intent(IN)  :: F, U
    real(real_kind), pointer, dimension(:)     :: U_Bound, XC_Bound

    ! local variables
    real(real_kind), dimension(ndim)  :: u1, u2, fj, Gvec
    real(real_kind)                   :: ndotudiff, normtrac
    integer(KIND = int_kind)          :: inode, nnum, yindex, gnum, gindex, &
                                         combination, gap_pos

    ! Initialize to an invalid value
    Node_Displacement_BC%Lambda(:,:) = -1.0

    ! Loop over all displacement BC nodes
    do inode = 1,SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       yindex = (nnum -1)*ndim
       u1 = U(yindex+1:yindex+ndim) + Vertex(nnum)%Coord
       fj = F(yindex+1:yindex+ndim)
       combination = Node_Displacement_BC%Combination(inode)
       select case(combination)
       ! At least one normal constraint 
       case(ONE_NORM_CONST, TWO_NORM_CONST, THREE_NORM_CONST, &
            ONE_D_ONE_NC, TWO_D_ONE_NC, ONE_D_TWO_NC)
          ! The gap normal vector and gap node number are in a different position, 
          ! depending  on the number of specified displacements
          select case(combination)
          case(ONE_NORM_CONST, TWO_NORM_CONST, THREE_NORM_CONST)
             gap_pos = 1
          case(ONE_D_ONE_NC, ONE_D_TWO_NC)
             gap_pos = 2
          case(TWO_D_ONE_NC)
             gap_pos = 3
          end select
          Gvec = Node_Displacement_BC%Normal(inode,gap_pos,:)
          gnum = Node_Displacement_BC%Gap_Node(inode,gap_pos)
          gindex = (abs(gnum) - 1) * ndim
          if (gnum > 0) then
             u2 = U(gindex+1:gindex+ndim) + Vertex(gnum)%Coord
          else
             u2 = U_Bound(gindex+1:gindex+ndim) + XC_Bound(gindex+1:gindex+ndim) 
          end if
          ndotudiff = SUM((u2-u1)*Gvec)
          normtrac = -SUM(fj*Gvec) * cscale(ndim*(nnum))/Node_Displacement_BC%Area(inode,gap_pos)
          Node_Displacement_BC%Gap_Disp(inode,gap_pos) = ndotudiff
          Node_Displacement_BC%Normal_Traction(inode,gap_pos) = normtrac

          if (Node_Displacement_BC%BC_Type(inode,gap_pos) == CONTACT) then
             Node_Displacement_BC%Lambda(inode,gap_pos) = GET_LAMBDA(ndotudiff, normtrac)
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == NORMAL_CONSTRAINT) then
             ! We have a normal constraint condition, not contact
             Node_Displacement_BC%Lambda(inode,gap_pos) = 1.0
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == FREE_INTERFACE) then
             ! We have a free interface
             Node_Displacement_BC%Lambda(inode,gap_pos) = 0.0
          end if
       end select
      ! At least two normal constraints
       select case(combination)
       case(TWO_NORM_CONST, THREE_NORM_CONST, ONE_D_TWO_NC)
          ! The gap normal vector is in a different position, depending 
          ! on the number of specified displacements
          select case(combination)
          case(TWO_NORM_CONST, THREE_NORM_CONST)
             gap_pos = 2
          case(ONE_D_TWO_NC)
             gap_pos = 3
          end select
          Gvec = Node_Displacement_BC%Normal(inode,gap_pos,:)
          gnum = Node_Displacement_BC%Gap_Node(inode,gap_pos)
          gindex = (abs(gnum) - 1) * ndim
          if (gnum > 0) then
             u2 = U(gindex+1:gindex+ndim) + Vertex(gnum)%Coord
          else
             u2 = U_Bound(gindex+1:gindex+ndim) + XC_Bound(gindex+1:gindex+ndim)
          end if
          ndotudiff = SUM((u2-u1)*Gvec)
          normtrac = -SUM(fj*Gvec) * cscale(ndim*(nnum))/Node_Displacement_BC%Area(inode,gap_pos)
          Node_Displacement_BC%Gap_Disp(inode,gap_pos) = ndotudiff 
          Node_Displacement_BC%Normal_Traction(inode,gap_pos) = normtrac
               
          if (Node_Displacement_BC%BC_Type(inode,gap_pos) == CONTACT) then
             Node_Displacement_BC%Lambda(inode,gap_pos) = GET_LAMBDA(ndotudiff, normtrac)
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == NORMAL_CONSTRAINT) then
             ! We have a normal constraint condition, not contact
             Node_Displacement_BC%Lambda(inode,gap_pos) = 1.0
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == FREE_INTERFACE) then
             ! We have a free interface
             Node_Displacement_BC%Lambda(inode,gap_pos) = 0.0
          end if
       end select

       ! Three normal constraints
       select case(combination)
       case(THREE_NORM_CONST)
          gap_pos = 3
          gnum = Node_Displacement_BC%Gap_Node(inode,gap_pos)
          Gvec = Node_Displacement_BC%Normal(inode,gap_pos,:)
          gindex = (abs(gnum) - 1) * ndim
          if (gnum > 0) then
             u2 = U(gindex+1:gindex+ndim) + Vertex(gnum)%Coord
          else
             u2 = U_Bound(gindex+1:gindex+ndim) + XC_Bound(gindex+1:gindex+ndim)
          end if
          ndotudiff = SUM((u2-u1)*Gvec)
          normtrac = -SUM(fj*Gvec) * cscale(ndim*(nnum))/Node_Displacement_BC%Area(inode,gap_pos)
          Node_Displacement_BC%Gap_Disp(inode,gap_pos) = ndotudiff
          Node_Displacement_BC%Normal_Traction(inode,gap_pos) = normtrac
               
          if (Node_Displacement_BC%BC_Type(inode,gap_pos) == CONTACT) then
             Node_Displacement_BC%Lambda(inode,gap_pos) = GET_LAMBDA(ndotudiff, normtrac)
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == NORMAL_CONSTRAINT) then
             ! We have a normal constraint condition, not contact
             Node_Displacement_BC%Lambda(inode,gap_pos) = 1.0
          else if (Node_Displacement_BC%BC_Type(inode,gap_pos) == FREE_INTERFACE) then
             ! We have a free interface
             Node_Displacement_BC%Lambda(inode,gap_pos) = 0.0
          end if
       end select
    end do
! This may be needed for contact functions using the force at a node
! The intent is to ensure that the lambda value is the same for contact
! node pairs.  This may be out of sync with the current code, and could be 
! expensive, especially if we try to get off-processor values.
! Uncomment Starting here...
!    ! Take A Second Pass To Force lambda To Be Symmetric
!    ! Assume That lambda Could Have Any Value Between 0 And 1.
!    do Inode = 1,Size(Node_Displacement_BC%node)
!       Nnum = Node_Displacement_BC%Node(Inode)
!       glambda = -1.0
!       do Idim = 1,Ndim
!          Gnum = Node_Displacement_BC%Gap_Node(Inode,Idim)
!          ! Skip Gap_node Elements Without Entries Or That Are Off-Processor
!          If(Gnum > 0) Then
!             Jnode = Nbc_index(Node_Displacement_BC%Gap_Node(Inode,Idim))
!             lambda =  Node_Displacement_BC%Lambda(Inode,Idim)
!             ! Find The Matching Node 
!             do Jdim = 1,Ndim
!                ! If The Node Matches...
!                If (Node_Displacement_BC%Gap_Node(Jnode,Jdim) == Nnum) Then
!                   Gvec = Node_Displacement_BC%Normal(Inode,Idim,:)
!                   Mgvec = Node_Displacement_BC%Normal(Jnode,Jdim,:)
!                   Ndot = Sum(Gvec * Mgvec)
!                   ! And The Vector Matches (Is The Opposite Sign)
!                   If ((Ndot < -0.9999) .And. (Ndot > -1.0001)) Then
!                      glambda = Node_Displacement_BC%Lambda(Jnode,Jdim)
!                      Exit
!                   End If
!                End If
!             End do 
!            If (glambda == -1.0) Then
!                Call Punt((/' Gap Node lambda Not Found'/), 'Contact_function')
!             End If
!             ! Consistency Check
!             If ((lambda < 0.0) .Or.(glambda < 0.0)) & 
!                    Call Punt((/' Invalid lambda '/), 'Contact_function')
!             ! Set Both lambdas To The Larger Value, Erring On The Side Of Forcing Contact
!             ! (If One Of The Values Is Zero And The Other Is Not, Then The Displacement 
!             ! Criterion Must Be In Favor Of Contact, And The Difference Is Due To A Force Mismatch.)
!             Node_Displacement_BC%Lambda(Inode,Idim) = Merge(lambda, glambda,(lambda >= glambda))!
!             Node_Displacement_BC%Lambda(Jnode,Jdim) = Node_Displacement_BC%Lambda(Inode,Idim)
!          End If
!       End do
!    End do

  END SUBROUTINE CONTACT_FUNCTION

  FUNCTION GET_LAMBDA(ndotudiff, normtrac) result (lam)

    real(real_kind), intent(IN) :: ndotudiff, normtrac
    real(real_kind)             :: lam

    ! local variables
    
! Simple polynomial of displacement only
!    if (ndotudiff <= 0.0) then
!       lam = 1.0
!    else if (ndotudiff >= contact_distance) then
!       lam = 0.0
!    else
!       lam = 2.0 * (ndotudiff/contact_distance - 1.0)**3 + 3.0 * (ndotudiff/contact_distance - 1.0)**2
!    end if

! Simple polynomial of displacement and normal traction
    ! Definitely in contact - (s < 0) and normal traction is compressive
    if ((ndotudiff <= 0.0) .and. (normtrac <=0.0)) then
       lam = 1.0
    ! Definitely not in contact
    else if ((ndotudiff >= contact_distance) .or. (normtrac >= contact_norm_trac)) then
       lam = 0.0
    ! s is small but positive, and normal traction is 0 or negative
    else if (normtrac <= 0.0) then
       lam = 2.0 * (ndotudiff/contact_distance - 1.0)**3 + 3.0 * (ndotudiff/contact_distance - 1.0)**2
    ! s is 0 or negative, and normal traction is small but tensile
    else if (ndotudiff <= 0.0) then
       lam = 2.0 * (normtrac/contact_norm_trac - 1.0)**3 + 3.0 * (normtrac/contact_norm_trac - 1.0)**2
    else
       lam = (2.0 * (ndotudiff/contact_distance - 1.0)**3 + 3.0 * (ndotudiff/contact_distance - 1.0)**2) * &
             (2.0 * (normtrac/contact_norm_trac - 1.0)**3 + 3.0 * (normtrac/contact_norm_trac - 1.0)**2)
    end if


  END FUNCTION GET_LAMBDA

  SUBROUTINE MECH_PRECOND_DISP_CONSTRAINTS(A_Elas)
    !=============================================================================
    !
    ! Apply nodal constraints in the preconditioning matrix such as interface and 
    ! contact constraints with a penalty method.
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================

    use mech_bc_data_module
    use parameter_module,     only: ndim
    use mesh_module,          only: Vertex_Ngbr_All
    implicit none

    ! Arguments
    type(real_var_vector), pointer, dimension(:) :: A_Elas

    ! Local variables
    integer                                      :: status
    Integer(KIND = int_kind)                     :: inode, yindex, nnum, idim, jdim, &
                                                    nmax, n, nindex
    real(real_kind),dimension(ndim)              :: Nvec1, Tvec1
    real(real_kind)                              :: atmp, n1dotf, t1dotf
    
    real(real_kind), pointer, dimension(:)       :: A_Vec
    real(real_kind), allocatable, dimension(:,:) :: A_Vec2, A_Save
  
    penalty = contact_penalty

    ! Apply projections of force and displacement vectors, making use of the 
    ! fact that the ith component of the vector [n n_T] u == n_i (n dot u)
    do inode = 1,SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       yindex = (nnum -1)*ndim
       nmax =   SIZES(Vertex_Ngbr_All(nnum))

       ! Working copy of the coefficients of the ndim equations to be modified
       allocate(A_Vec2(ndim,(nmax+1)*ndim), stat = status)
       if (status /= 0) call TLS_panic ('MECH_PRECOND_DISP_CONSTRAINTS: allocation error: A_Vec2')
       ! Also keep a copy of the unmodified coefficients
       allocate(A_Save(ndim,(nmax+1)*ndim), stat = status)
       if (status /= 0) call TLS_panic ('MECH_PRECOND_DISP_CONSTRAINTS: allocation error: A_Save')
       do idim = 1,ndim
          A_Vec => FLATTEN(A_Elas(yindex + idim))
          A_Vec2(idim,:) = A_Vec
       end do

       ! Rearrange the first three components to be x,y,z, instead of diagonal in the first 
       ! element in the array
       atmp = A_Vec2(2,2)
       A_Vec2(2,2) = A_Vec2(2,1)
       A_Vec2(2,1) = atmp
       atmp = A_Vec2(3,1)
       A_Vec2(3,1) = A_Vec2(3,3)
       A_Vec2(3,3) = atmp

       A_Save = A_Vec2

       select case(Node_Displacement_BC%Combination(inode))
       case(ONE_DISPLACEMENT, ONE_D_ONE_NC, ONE_D_TWO_NC)
          ! Normal vector
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Accumulate nn^T(f_j) - force term
          do n = 1,nmax + 1
             nindex = (n-1) * ndim
             do jdim = 1,ndim
                n1dotf = 0.0
                n1dotf = SUM(Nvec1(:) * A_Vec2(:,nindex + jdim))
                do idim = 1,ndim
                   A_Vec2(idim,nindex + jdim) = A_Vec2(idim,nindex + jdim) - n1dotf * Nvec1(idim)
                end do
             end do
          end do
          ! Specify normal displacement for each DOF for this node (inode)
          do idim = 1,ndim
             ! Penalize the displacement normal to the interface (need all displacement
             ! components.
             do jdim = 1,ndim
                A_Vec2(idim,jdim) = A_Vec2(idim,jdim) - penalty * Nvec1(jdim) * Nvec1(idim)
             end do
          end do
       case(TWO_DISPLACEMENTS, TWO_D_ONE_NC)
          ! Tangent vector t \propto n1 x n2
          Tvec1 = Node_Displacement_BC%Vector(inode,2,:)
          ! Accumulate t dot f
          do n = 1,nmax + 1
             nindex = (n-1) * ndim
             ! Tangential component of f
             do jdim = 1,ndim
                t1dotf = 0.0
                t1dotf = SUM(Tvec1(:) * A_Vec2(:,nindex + jdim))
                do idim = 1,ndim
                   A_Vec2(idim,nindex + jdim) = t1dotf * Tvec1(idim)
                end do
             end do
          end do
          ! Add in-plane component of u
          do idim = 1,ndim
             do jdim = 1,ndim
                A_Vec2(idim,jdim) = A_Vec2(idim,jdim) + penalty * Tvec1(jdim) * Tvec1(idim)
                if (idim == jdim) A_Vec2(idim,jdim) = A_Vec2(idim,jdim) - penalty
             end do
          end do
       case(THREE_DISPLACEMENTS)
          ! We already have a unique displacement vector in global coordinates
          ! The displacement values are set in the right hand side
          A_Vec2 = 0.0
          A_Vec2(1,1) = - penalty 
          A_Vec2(2,2) = - penalty 
          A_Vec2(3,3) = - penalty 

       end select

       ! Restore diagonals in element 1
       atmp = A_Vec2(2,2)
       A_Vec2(2,2) = A_Vec2(2,1)
       A_Vec2(2,1) = atmp
       atmp = A_Vec2(3,1)
       A_Vec2(3,1) = A_Vec2(3,3)
       A_Vec2(3,3) = atmp
       ! Copy A_Vec2 back to the preconditioner matrix
       do idim = 1,ndim
          A_Vec => FLATTEN(A_Elas(yindex + idim))
          A_Vec = A_Vec2(idim,:)
       end do
       Nullify(A_Vec)
       deallocate(A_Vec2)
       deallocate(A_Save)
    end do

  END SUBROUTINE MECH_PRECOND_DISP_CONSTRAINTS

  !! Copied from MeshSupport
  pure function cross_product (a, b) result (axb)
    use kind_module,          only: real_kind
    real(real_kind), intent(in) :: a(:), b(:)
    real(real_kind)             :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  SUBROUTINE FACE_GAP_INITIALIZE()
    !=============================================================================
    !
    ! Set up data for using  gap opening displacements at nodes to calculate a 
    ! gap displacement value for the cell faces.
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================

    use kind_module,          only: int_kind, real_kind
    use mech_bc_data_module
    use bc_operations
    use parameter_module,     only: nvf
    use pgslib_module,        only: PGSLib_GLOBAL_ANY
    implicit none

    TYPE(BC_Operator), POINTER                 :: HTC_GAP_Operator
    TYPE(BC_Atlas),    POINTER                 :: HTC_GAP_Atlas

    INTEGER(int_kind), POINTER, DIMENSION(:)   :: FaceList
    INTEGER(int_kind), POINTER, DIMENSION(:)   :: CellList
    INTEGER(int_kind), POINTER, DIMENSION(:)   :: OffsetList
    REAL(real_kind), POINTER, DIMENSION(:,:)   :: ValueList

    Integer(KIND = int_kind)                   :: numfaces, iface, jnode, l, &
                                                  finterface, status
    Integer(KIND = int_kind),dimension(nvf)    :: v

    numfaces = 0

    ! Get the HTC_GAP operator
    HTC_GAP_Operator => BC_Spec_Get_Operator(TEMPERATURE_BC, BC_HTC_GAP_Op)

    ! We want to traverse the whole atlas.
    HTC_GAP_Atlas => BC_OP_Get_Atlas(HTC_GAP_Operator)

    numfaces = DATA_SIZE(HTC_GAP_Atlas)
    IF (PGSLib_Global_ANY(numfaces > 0)) THEN

       ! We are using whole arrays, so faster to point at them than get a private copy.
       FaceList     => BC_Get_Face(HTC_GAP_Atlas)
       CellList     => BC_Get_Cell(HTC_GAP_Atlas)
       OffsetList   => BC_Get_Offset(HTC_GAP_Atlas)
       ValueList    => BC_Get_Values(HTC_GAP_Atlas)

       allocate(Face_Gap(numfaces),stat = status)
       if (status /= 0) call TLS_panic ('FACE_GAP_INITIALIZE: allocation error: Face_Gap')
       allocate(Face_Node(2,nvf,numfaces),stat = status)
       if (status /= 0) call TLS_panic ('FACE_GAP_INITIALIZE: allocation error: Face_Node')
       Face_Gap = 0.0
       Face_Node = 0

       do iface = 1,numfaces
          select case (FaceList(iface))
          case (1)
             v(1)=3; v(2)=4; v(3)=8; v(4)=7
          case (2)
             v(1)=1; v(2)=2; v(3)=6; v(4)=5
          case (3)
             v(1)=4; v(2)=1; v(3)=5; v(4)=8
          case (4)
             v(1)=2; v(2)=3; v(3)=7; v(4)=6
          case (5)
             v(1)=1; v(2)=2; v(3)=3; v(4)=4
          case (6)
             v(1)=7; v(2)=8; v(3)=5; v(4)=6
          end select
          finterface = Interface_ID(FaceList(iface), CellList(iface))
          if (finterface /= 0) then
             do l = 1,SIZE(Interface_List)
                if (Interface_List(l) == finterface) &
                     Face_Node(2,:,iface) = l
             end do
          end if

          do jnode = 1, nvf
             Face_Node(1,jnode,iface) = v(jnode)
          end do
          ! Move value data from (1:3) to (2:4) so that we can use the first 
          ! element for the current heat transfer coefficient.  The first 
          ! element is set to zero, but will be changed by the next call.
          ValueList(4,OffsetList(iface)) = ValueList(3,OffsetList(iface)) 
          ValueList(3,OffsetList(iface)) = ValueList(2,OffsetList(iface)) 
          ValueList(2,OffsetList(iface)) = ValueList(1,OffsetList(iface)) 
          ValueList(1,OffsetList(iface)) = 0.0
       end do
    END IF
       
  END SUBROUTINE FACE_GAP_INITIALIZE

  SUBROUTINE FACE_GAP_UPDATE()
    !=============================================================================
    !
    ! Use gap opening displacements at nodes to calculate a gap displacement value 
    ! for the cell faces.
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================

    use kind_module,          only: int_kind, real_kind
    use bc_operations
    use mech_bc_data_module
    use parameter_module,     only: nvf, nvc, ncells
    use pgslib_module,        only: PGSLib_GLOBAL_ANY
    use gs_module,            only: EN_Gather
    implicit none

    TYPE(BC_Operator), POINTER                 :: HTC_GAP_Operator
    TYPE(BC_Atlas),    POINTER                 :: HTC_GAP_Atlas

    Integer(KIND = int_kind)                   :: nv, iface, numfaces, iint, &
                                                  status, vrtx, icell
    real(kind=real_kind),dimension(nvf)        :: g = 0.0

    INTEGER(int_kind), POINTER, DIMENSION(:)   :: CellList
    INTEGER(int_kind), POINTER, DIMENSION(:)   :: OffsetList
    REAL(real_kind), POINTER, DIMENSION(:,:)   :: ValueList

    real(real_kind), pointer, dimension(:,:)     :: ENtemp
    real(real_kind), pointer, dimension(:,:,:)   :: EN_Gap
     
    NULLIFY(ENtemp)
    NULLIFY(EN_Gap)

    ! Get the HTC_GAP operator
    HTC_GAP_Operator => BC_Spec_Get_Operator(TEMPERATURE_BC, BC_HTC_GAP_Op)

    ! We want to traverse the whole atlas.
    HTC_GAP_Atlas => BC_OP_Get_Atlas(HTC_GAP_Operator)

    numfaces = DATA_SIZE(HTC_GAP_Atlas)
    IF (PGSLib_Global_ANY(numfaces > 0)) THEN

       ! Gather gap displacements of off-processor neighboring nodes
       allocate(ENtemp(nvc, ncells),stat = status)
       if (status /= 0) call TLS_panic ('FACE_GAP_UPDATE: allocation error: ENtemp')    
       allocate(EN_Gap(nvc, ncells, SIZE(Node_Gap,2)),stat = status)
       if (status /= 0) call TLS_panic ('FACE_GAP_UPDATE: allocation error: EN_Gap')    
       do iint = 1,SIZE(Node_Gap,2)
          call EN_GATHER(ENtemp, Node_Gap(:,iint))
          EN_Gap(:,:,iint) = ENtemp(:,:)
       end do

       ! We are going to change the values
       CellList     => BC_Get_Cell(HTC_GAP_Atlas)
       OffsetList   => BC_Get_Offset(HTC_GAP_Atlas)
       ValueList    => BC_Get_Values(HTC_GAP_Atlas)

       do iface = 1,SIZE(Face_Gap)
          do nv = 1, nvf
             if (Face_Node(2,nv,iface) /= 0) then
                icell = CellList(iface)
                vrtx = Face_Node(1,nv,iface)
                g(nv) = EN_Gap(vrtx,icell,Face_Node(2,nv,iface))
             else
                g(nv) = 0.0
             end if
          end do
          Face_Gap(iface) = SUM(g(:))/nvf

          if (Face_Gap(iface) > ValueList(4,OffsetList(iface))) then
             ValueList(1,OffsetList(iface)) = ValueList(3,OffsetList(iface))
          else
             ValueList(1,OffsetList(iface)) = ValueList(2,OffsetList(iface))
          end if
       end do

       deallocate(ENtemp)
       deallocate(EN_Gap)
       
    end IF

  END SUBROUTINE FACE_GAP_UPDATE

end Module SOLID_MECH_CONSTRAINTS

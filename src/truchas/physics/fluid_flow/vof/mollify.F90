MODULE MOLLIFY
  !=======================================================================
  ! Purpose:
  !    Set up derived types and fill with the cell number of neighbors
  !    within each convolution
  !
  ! Author(s): M. Williams (MST-8, mww@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: MOLLIFY_CONV_SAVEMEM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  type MOLLIFY_TYPE
 
     ! Cell numbers of cells in mollification domain
     integer, DIMENSION(:), pointer :: ngbr_id

     ! Distance to cell centroids of cells in mollification domain
     real(r8), DIMENSION(:,:), pointer :: ngbr_dist_vector

     ! Cell volumes of cells in mollification domain
     real(r8), DIMENSION(:), pointer :: ngbr_vol
 
  end type MOLLIFY_TYPE

  ! Variables in NUMERICS namelist
  real(r8), public, save :: interface_smoothing_length

CONTAINS

  SUBROUTINE MOLLIFY_CONV_SAVEMEM (Scal_Field,Moll_Scal,N_x,N_y,N_z, matl_no)
    !=======================================================================
    !  Purpose:
    !    Returns the mollified scalar field 
    !    This routine is typically used to convolve volume fraction, though
    !    it is meant to be generic as to convolve any given scalar field.
    !    The `convolved normals' algorithm is included here so we don't have
    !    to allocate and deallocate all of the ncells_tot length arrays.
    !=======================================================================
       use cutoffs_module,          only: cutvof
       use mesh_module,             only: Mesh,Cell, degenerate_face
       use pgslib_module,           only: PGSLIB_BCast, PGSLIB_Collate &
                                      ,PGSLib_GLOBAL_Maxval
       use parameter_module,        only: ncells, ncells_tot, ndim, nfc, nmat
       use var_vector_module

       real(r8), dimension(ncells), INTENT(IN) :: Scal_Field
       real(r8), dimension(ncells), INTENT(INOUT) :: Moll_Scal
       real(r8), dimension(nmat,ncells), INTENT(INOUT) :: N_x,N_y,N_z
       integer, INTENT(IN) :: matl_no
 
!... declare local variables
      integer :: I, I0, K, C, dim_cnt, memstat
      real(r8) :: dist_x, dist_y, dist_z, KF, KdxF, KDyF, KdzF
      real(r8), dimension(ndim) :: cent, dist
      logical, allocatable, dimension(:) :: Hit_Good, Mask
      logical, save :: first_time = .true.
      integer,  pointer, save, dimension(:,:)   :: All_Ngbrs_All, All_Ngbr
      real(r8), pointer, save, dimension(:)     :: Tot_Vol, Tot_Scal
      real(r8), pointer, save, dimension(:,:)   :: All_centroid
      real(r8), pointer, save, dimension(:,:,:) :: Tot_Face_Centroid
      ! variable to find smooth interface indicator function
      real(r8) :: proc_integ, tot_integ, DEL_INT

     ALLOCATE(Hit_Good(ncells_tot),Mask(ncells),STAT = memstat)
     call TLS_fatal_if_any ((memstat /= 0), 'MOLLIFY_CONV_SAVEMEM: Memory allocation error for logical arrays')

! ...
! ... Find maximum Ngbr_Cells_All, used to allocate rectangular all_ngbr array
! ... 

 if(first_time)then
       ALLOCATE(All_Ngbrs_All(ncells_tot,nfc),All_Ngbr(ncells,nfc),   &
       Tot_Vol(ncells_tot),Tot_Scal(ncells_tot),                      &
       Tot_Face_Centroid(ncells_tot,ndim,nfc),                        &
       All_centroid(ncells_tot,ndim), STAT = memstat)
       call TLS_fatal_if_any ((memstat /= 0), 'MOLLIFY_CONV_SAVEMEM: Memory allocation error')

! ... Initialize arrays
       All_Ngbr = 0

! ...
! ... We need Ngbr_Cells_All of ALL cells in the domain
! ... fist pack on processor Ngbr_Cells_All into All_Ngbr, then collate & dist
! ... We'll build the ragged array later, right now punt

     Do I=1,ncells
      Do k=1,nfc
       All_Ngbr(I,k) = Mesh(I)%ngbr_Cell_orig(k)
      End do
     End do


    Do I=1,nfc
     call PGSLIB_Collate(All_Ngbrs_All(:,i),All_Ngbr(:,i))
    End do
    Do I=1,nfc
     call PGSLIB_Bcast(All_Ngbrs_All(:,i))
    Enddo

! ...
! ... Now to Collate and Broadcast All Cell Centroids,Cell Volumes,Scalar Field
    call PGSLIB_Collate(Tot_Vol(:), Cell(:)%Volume)
    call PGSLIB_Bcast(Tot_Vol)
     do dim_cnt = 1,ndim
         call PGSLIB_Collate(All_Centroid(:,dim_cnt),Cell(:)%Centroid(dim_cnt))
     end do
     do dim_cnt = 1,ndim
         call PGSLIB_Bcast(All_Centroid(:,dim_cnt))
     end do

! ...
! ... Now to Collate and Broadcast All Face Centroids

    do k=1,nfc
     do dim_cnt = 1,ndim
         call PGSLIB_Collate(Tot_Face_Centroid(:,dim_cnt,k),Cell(:)%Face_Centroid(dim_cnt,k))
     end do
    end do
    do k=1,nfc
     do dim_cnt = 1,ndim
         call PGSLIB_Bcast(Tot_Face_Centroid(:,dim_cnt,k))
     end do
    enddo
    first_time = .false.
end if
! ...
! ... Now to Collate and Broadcast Scalar Field (volume fractions)
    call PGSLIB_Collate(Tot_Scal(:), Scal_Field(:))
    call PGSLIB_Bcast(Tot_Scal)

! ...
! ... SWEEP finds the value of the convolution in each cell by using a
! .. recursive subroutine which examines cell neighbors and then the
     ! ... neighbors of neighbors
! ...
      del_int = 0.8* interface_smoothing_length
      proc_integ = 0.0

    SWEEP_VOF: Do I0=1,ncells

! ... find the value of the numeric integration centered at the centroid 
! ... of the unit cube  domain 
      proc_integ = proc_integ +               &
              Kern(0.0-Cell(I0)%Centroid(1), 0.0-Cell(I0)%Centroid(2),    &
                   0.0-Cell(I0)%Centroid(3), del_int)*Cell(I0)%Volume
!              Kern(0.525-Cell(I0)%Centroid(1), 0.525-Cell(I0)%Centroid(2),    &
!                   0.525-Cell(I0)%Centroid(3), del_int)*Cell(I0)%Volume
! ... Initialize Hit array and convolution value
! ... and set the centroid for the convolution
     Hit_Good = .false.
     KF = 0.0

      Do dim_cnt = 1,ndim
         cent(dim_cnt) =  Cell(I0)%Centroid(dim_cnt)
      End Do

! .. begin recursive loop

   RECURS_VOF:  DO I = 1,nfc
           C = Mesh(I0)%Ngbr_Cell_orig(I)
            IF(C == DEGENERATE_FACE) cycle RECURS_VOF
            IF(C == 0)then
! ...
! ... NOTE: the BC's are only for homogenous neumann on the unit cube
! ...
             IF(I == 1)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = -0.005 + (-0.005 -cent(dim_cnt)) - cent(dim_cnt)
               Case (2)
                dist_y =  cent(dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z =  cent(dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End IF
             IF(I == 2)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = 0.005 + (0.005- cent(dim_cnt)) - cent(dim_cnt)
               Case (2)
                dist_y = cent(dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = cent(dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF(I == 3)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = cent(dim_cnt) - cent(dim_cnt)
               Case (2)
!                dist_y =  - cent(dim_cnt) - cent(dim_cnt)
                dist_y =  -0.00025 + (-0.00025 - cent(dim_cnt)) - cent(dim_cnt)
               Case (3)
                dist_z = cent(dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF(I == 4)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = cent(dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = 0.00025 + (0.00025 - cent(dim_cnt)) - cent(dim_cnt)
               Case (3)
                dist_z = cent(dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End if
             If( I == 5)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt) 
               Case (1)
                dist_x =  cent(dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y =  cent(dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = -0.005 + (-0.005 - cent(dim_cnt)) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF( I == 6)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = cent(dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = cent(dim_cnt) - cent(dim_cnt)
               Case (3)
                 dist_z = 0.005 + (0.005 - cent(dim_cnt)) - cent(dim_cnt)
               End Select
              End Do
             EndIf
               KF = KF + Scal_Field(I0)*                                  &
                  Kern(dist_x,dist_y, dist_z, del_int)*Cell(I0)%Volume

            ELSE
            IF(C > 0 .AND. .not. Hit_Good(c))then
             Hit_Good(c) = .true.
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_centroid(c,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = all_centroid(c,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_centroid(c,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             IF(ndim == 2)dist_z = 0.0
              IF(SQRT(dist_x**2+dist_y**2+dist_z**2) < del_int)then
               KF = KF + Tot_Scal(c)*                                            &
                  Kern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)

               CALL RECURS_CONV(C,Tot_Vol,Tot_Scal,All_Ngbrs_All, All_centroid,  &
                                    cent,Hit_Good,KF, del_int)
              END IF
            END IF
            END IF
          ENDDO RECURS_VOF

      Moll_Scal(I0) = KF
    End Do SWEEP_VOF

!      tot_integ = PGSLIB_Global_Sum(proc_integ)
      tot_integ = PGSLIB_Global_MaxVal(Moll_Scal(:))

! ... Now to find the unit normals to the interface

! ...
! ... Set up a Mask to find normals - for interface or off-interface normals
! ... 
!      IF(surface_tension_model == 'convolution model')then
!        Mask = (Moll_Scal(:)  < tot_integ-0.001 .AND. &
!                     Moll_Scal(:)  > .001)
!         Mask = (Moll_Scal(:)  < 0.9995*tot_integ .AND. &
!                     Moll_Scal(:)  > 0.0005*tot_integ)
!      ELSE
         Mask = (Scal_Field < 1.0-cutvof .AND. Scal_Field > cutvof)
!      END IF
! ...
! ... Now sweep to find interface normals where Mask=.true.
! ... find normals everywhere - do the masking in the surf10 algorithm
      Mask = .true.
     DEL_INT = 1.2 *interface_smoothing_length
    SWEEP_NORM: Do I0=1,ncells
    IF(MASK(I0))then

! ... Initialize Hit array and convolution value
! ... and set the centroid for the convolution
     Hit_Good = .false.
     KDxF = Scal_Field(I0)*DzKern(0.D0,0.D0,0.D0,del_int)*Cell(I0)%Volume
     KDyF = Scal_Field(I0)*DzKern(0.D0,0.D0,0.D0,del_int)*Cell(I0)%Volume
     KDzF = Scal_Field(I0)*DzKern(0.D0,0.D0,0.D0,del_int)*Cell(I0)%Volume
     Hit_Good(I0) = .true.

      Do dim_cnt = 1,ndim
         cent(dim_cnt) =  Cell(I0)%Centroid(dim_cnt)
      End Do

! .. begin recursive loop

     RECURS_NORM_LOOP: DO I = 1,SIZE(Mesh(I0)%Ngbr_Cell)
           C = Mesh(I0)%Ngbr_Cell_orig(I)
            IF(C == DEGENERATE_FACE) cycle RECURS_NORM_LOOP
            IF(C == 0)then
! ...
! ... NOTE: the BC's are only for homogenous neumann on the unit cube
! ...
              dist(:) = 2.D0*(Cell(I0)%Face_Centroid(:,I) - cent(:))

              dist_x = dist(1)
              dist_y = dist(2)
              dist_z = dist(3)

              KDxF = KDxF + Scal_Field(I0)*                                   &
                  DxKern(dist_x,dist_y, dist_z, del_int)*Cell(I0)%Volume
              KDyF = KDyF + Scal_Field(I0)*                                   &
                  DyKern(dist_x,dist_y, dist_z, del_int)*Cell(I0)%Volume
              KDzF = KDzF + Scal_Field(I0)*                                   &
                  DzKern(dist_x,dist_y, dist_z, del_int)*Cell(I0)%Volume
            ELSE
            IF(C > 0 .AND. .not. Hit_Good(c))then
             Hit_Good(c) = .true.
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_centroid(c,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = all_centroid(c,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_centroid(c,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             IF(ndim == 2)dist_z = 0.0

               KDxF = KDxF + Tot_Scal(c)*                                      &
                  DxKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
               KDyF = KDyF + Tot_Scal(c)*                                      &
                  DyKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
               KDzF = KDzF + Tot_Scal(c)*                                      &
                  DzKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)

             CALL RECURS_NORM(C,Tot_Vol,Tot_Scal,All_Ngbrs_All, All_centroid,  &
                     Tot_Face_Centroid, cent,Hit_Good,KDxF,KDyF,KdzF, del_int)
              END IF
            END IF
          ENDDO RECURS_NORM_LOOP

      N_x(matl_no,I0) = KDxF
      N_y(matl_no,I0) = KDyF
      N_z(matl_no,I0) = KDzF
! ... use dist_x as a temp to normalize the normal length
      dist_x = SQRT(N_x(matl_no,I0)**2+N_y(matl_no,I0)**2+N_z(matl_no,I0)**2)
       IF(dist_x > 1.D-09)then
         N_x(matl_no,I0) = N_x(matl_no,I0)/dist_x
         N_y(matl_no,I0) = N_y(matl_no,I0)/dist_x
         N_z(matl_no,I0) = N_z(matl_no,I0)/dist_x
       END IF
    ELSE
! ... This ELSE condition is for Mask=.false.
      N_x(matl_no,I0) = 0.0
      N_y(matl_no,I0) = 0.0
      N_z(matl_no,I0) = 0.0

    ENDIF
    End Do SWEEP_NORM
    if(Associated(All_Ngbr))DEALLOCATE(All_Ngbr)
    DEALLOCATE(Hit_Good, Mask)

  END SUBROUTINE MOLLIFY_CONV_SAVEMEM

  RECURSIVE SUBROUTINE RECURS_CONV(Cell_no,Tot_vol, Tot_scal,all_ngbr,    &
                                             all_cent, cent,Hit_Array,KF,del_int)
    !=======================================================================
    !  Purpose:
    !    Recursive subroutine to calculate the value of the convolution of 
    !    Tot_Scal with a kernel centered at some given cell -
    !    This program is called from the CONV_MOLLIFY_SAMEMEM routine 
    !=======================================================================
      use mesh_module, only: DEGENERATE_FACE
      use parameter_module, only: ndim

      integer, INTENT(IN) :: Cell_no
      real(r8), dimension(:), INTENT(IN) :: Tot_vol
      real(r8), dimension(:), INTENT(IN) :: Tot_Scal
      integer,  dimension(:,:), INTENT(IN) :: all_ngbr
      real(r8), dimension(:,:), INTENT(IN) :: all_cent
      real(r8), dimension(:), INTENT(IN) :: cent
      logical, dimension(:), INTENT(INOUT) :: Hit_Array
      real(r8), INTENT(INOUT) :: KF
      real(r8), INTENT(IN) :: del_int

! ... Local variables
      integer :: J,C, dim_cnt
      real(r8) :: radius,dist_x,dist_y,dist_z

! ... 
! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are
! ... 

        RECURS_CONV_LOOP: Do J = 1,SIZE(all_ngbr(Cell_no,:))
         c = all_ngbr(Cell_no,J)
            IF(C == DEGENERATE_FACE) cycle RECURS_CONV_LOOP
            IF(C == 0)then
             IF(J == 1)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = -all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y =  all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z =  all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End IF
             IF(J == 2)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = 1.0 + (1.0- all_cent(Cell_no,dim_cnt)) - cent(dim_cnt)
               Case (2)
                dist_y = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF(J == 3)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y =  - all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF(J ==4)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = 1.0 + (1.0 - all_cent(Cell_no,dim_cnt)) - cent(dim_cnt)
               Case (3)
                dist_z = all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End if
             If(J == 5)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x =  all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y =  all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = - all_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
             End If
             IF(J == 6)then
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x =  All_cent(Cell_no,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = all_cent(cell_no,dim_cnt)  - cent(dim_cnt)
               Case (3)
                dist_z = 1.0 + (1.0 - all_cent(Cell_no,dim_cnt)) - cent(dim_cnt)
               End Select
              End Do
             EndIf
               KF = KF + Tot_Scal(cell_no)*                             &
                  Kern(dist_x,dist_y, dist_z, del_int)*Tot_vol(Cell_no)

            ELSE

          IF(c > 0 .AND. .not. Hit_Array(c))then
           Hit_Array(c) = .true.
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_cent(c,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = all_cent(c,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_cent(c,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
           if(ndim == 2)dist_z = 0.0
           radius = SQRT(dist_x**2+dist_y**2+dist_z**2)

            IF(radius < interface_smoothing_length)then
             KF = KF + Tot_Scal(c)*                                            &
                Kern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
             CALL RECURS_CONV(c,Tot_Vol,Tot_Scal,all_ngbr,all_cent,cent,       &
                              Hit_Array,KF, del_int)
            ENDIF
          ENDIF
          END IF
        End Do RECURS_CONV_LOOP

  END SUBROUTINE RECURS_CONV

  RECURSIVE SUBROUTINE RECURS_NORM(Cell_no,Tot_vol, Tot_scal,all_ngbr,all_cent  &
                        ,all_face_cent, cent,Hit_Array,KDxF,KDyF,KDzF, del_int)
    !=======================================================================
    !  Purpose:
    !    Recursive subroutine to calculate the unit normals to an interface
    !    by convolving Tot_Scal with the derivatives of a kernel centered at 
    !    some given cell -
    !    This program is called from the CONV_MOLLIFY_SAMEMEM routine 
    !=======================================================================
      use mesh_module, only: DEGENERATE_FACE
      use parameter_module, only: ndim
      integer, INTENT(IN) :: Cell_no
      real(r8), dimension(:), INTENT(IN) :: Tot_vol
      real(r8), dimension(:), INTENT(IN) :: Tot_Scal
      integer, dimension(:,:), INTENT(IN) :: all_ngbr
      real(r8), dimension(:,:), INTENT(IN) :: all_cent
      real(r8), dimension(:,:,:), INTENT(IN) :: all_face_cent
      real(r8), dimension(:), INTENT(IN) :: cent
      logical, dimension(:), INTENT(INOUT) :: Hit_Array
      real(r8), INTENT(INOUT) :: KDxF,KDyF,KDzF
      real(r8), INTENT(IN) :: DEL_INT
! ... Local variables
      integer :: J,C, dim_cnt
      real(r8) :: radius,dist_x,dist_y,dist_z
      real(r8), dimension(ndim) :: dist

! ... 
! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are
! ... 
       RECURS_NORM_LOOP: Do J = 1,SIZE(all_ngbr(Cell_no,:))
         c = all_ngbr(Cell_no,J)
            IF(C == DEGENERATE_FACE) cycle RECURS_NORM_LOOP
            IF(C == 0)then

             dist(:) = all_cent(Cell_no,:)-cent(:) + &
               2.D0*(All_Face_Cent(Cell_no,:,J)-all_cent(Cell_no,:))

              dist_x = dist(1)
              dist_y = dist(2)
              dist_z = dist(3)

            KDxF = KDxF + Tot_Scal(cell_no)*                                         &
                DxKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(cell_no)
            KDyF = KDyF + Tot_Scal(cell_no)*                                         &
                DyKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(cell_no)
            KDzF = KDzF + Tot_Scal(cell_no)*                                         &
                DzKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(cell_no)

            ELSE
          IF(c > 0 .AND. .not. Hit_Array(c))then
            Hit_Array(c) = .true.
              Do dim_cnt = 1,ndim
               Select Case (dim_cnt)
               Case (1)
                dist_x = all_cent(c,dim_cnt) - cent(dim_cnt)
               Case (2)
                dist_y = all_cent(c,dim_cnt) - cent(dim_cnt)
               Case (3)
                dist_z = all_cent(c,dim_cnt) - cent(dim_cnt)
               End Select
              End Do
           if(ndim == 2)dist_z = 0.0
           radius = SQRT(dist_x**2+dist_y**2+dist_z**2)

           IF(radius < del_int)then

            KDxF = KDxF + Tot_Scal(c)*                                         &
                DxKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
            KDyF = KDyF + Tot_Scal(c)*                                         &
                DyKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
            KDzF = KDzF + Tot_Scal(c)*                                         &
                DzKern(dist_x,dist_y, dist_z, del_int)*Tot_vol(c)
            CALL RECURS_NORM(c,Tot_Vol,Tot_Scal,all_ngbr,all_cent,all_face_cent,cent, &
                    Hit_Array,KDxF,KDyF,KDzF, del_int)
           ENDIF
           ENDIF
          ENDIF
        End Do RECURS_NORM_LOOP

  END SUBROUTINE RECURS_NORM

!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>
! Convolution Kernel and it's derivatives listed below
!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>

  real FUNCTION Kern(X,Y,Z,DELTA)
    use constants_module,  only: pi
    use parameter_module,  only: ndim
    real(r8) :: X,Y,Z,DELTA,A
    IF (ndim == 2) then
      A = 4.0/(DELTA**8*pi)
    ELSE
      A = 315.0/(64.0*DELTA**9*pi)
    END IF
    IF ((X)**2+(Y)**2+(Z)**2 <= DELTA**2) THEN
      Kern = A*(DELTA**2-((X)**2+(Y)**2+(Z)**2))**3
    ELSE
      Kern = 0.0
    END IF
  END FUNCTION Kern
      
  real FUNCTION DxKern(X,Y,Z,DELTA)
    use constants_module, only: pi
    use parameter_module, only: ndim
    real(r8) :: DELTA,X,Y,Z,A
    IF (ndim == 2) then
      A = 4.0/(DELTA**8*pi)
    ELSE
      A = 315.0/(64.0*DELTA**9*pi)
    END IF
    IF (X**2+Y**2+Z**2 <= DELTA**2) THEN
      DxKern = -6*A*X*(DELTA**2-X**2-Y**2-Z**2)**2
    ELSE
      DxKern = 0.0
    END IF 
  END FUNCTION DxKern
      
  real FUNCTION DyKern(X,Y,Z,DELTA)
    use constants_module, only: pi
    use parameter_module, only: ndim
    real(r8) :: X,Y,Z,DELTA,A
    IF (ndim == 2) then
      A = 4.0/(DELTA**8*pi)
    ELSE
      A = 315.0/(64.0*DELTA**9*pi)
    END IF
    IF (X**2+Y**2+Z**2 <= DELTA**2) THEN
      DyKern = -6*A*Y*(DELTA**2-X**2-Y**2-Z**2)**2
    ELSE
      DyKern = 0.0
    END IF
  END FUNCTION DyKern
      
  real FUNCTION DzKern(X,Y,Z,DELTA)
    use constants_module, only: pi
    use parameter_module, only: ndim
    real(r8) :: DELTA,X,Y,Z,A
    IF (ndim == 2) then
      A = 4.0/(DELTA**8*pi)
    ELSE
      A = 315.0/(64.0*DELTA**9*pi)
    END IF
    IF (X**2+Y**2+Z**2 <= DELTA**2) THEN
      DzKern = -6*A*Z*(DELTA**2-X**2-Y**2-Z**2)**2
    ELSE
      DzKern = 0.0
    END IF
  END FUNCTION DzKern

END MODULE MOLLIFY

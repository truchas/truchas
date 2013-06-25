MODULE KERNEL_INTERPOLATION_MODULE 
  !=======================================================================
  ! Purpose:
  !   Convolution and interpolation using the compact kernel 
  !   given in Rudman IJNMF 1998
  !  
  !   based on Matthew Williams mollify subroutines
  !
  !   Marianne Francois,CCS-2,LANL, mmfran@lanl.gov  
  ! 
  !=======================================================================
  use truchas_logging_services
  implicit none
 
  ! Private Module
  private

  ! Public Procedures and Types
  public :: KERN_CONVOLUTION_CENTER, KERN_INTERPOLATION_FACE 

  ! File Identifier
  !Character (LEN=*), Parameter :: version = &

CONTAINS
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE KERN_CONVOLUTION_CENTER(Scal_Center,Scal_Convol,d)

    !=======================================================================
    !  Purpose:
    !    Returns the convolved value of Scal_Center
    !    using Rudman 1998 kernel 
    !=======================================================================
    use constants_module,        only: zero 
    use kind_module,             only: real_kind, int_kind, log_kind
    use mesh_module,             only: Mesh,Cell, degenerate_face
    use pgslib_module,           only: PGSLIB_BCast, PGSLIB_Collate
    use parameter_module,        only: ncells, ncells_tot, ndim, nfc
    use var_vector_module

    IMPLICIT NONE

       real(real_kind), dimension(ncells),      INTENT(IN)    :: Scal_Center
       real(real_kind), dimension(ncells),      INTENT(INOUT) :: Scal_Convol 
       real(real_kind)                                        :: d
! local variables
      integer(int_kind)                    :: n,f,m,icell, memstat
      real(real_kind)                      :: dist_x, dist_y, dist_z,Ke,KF_N,KF_D 
      real(real_kind),dimension(ndim)      :: cent

      real(real_kind),allocatable,   dimension(:)     :: Tot_Vol, Tot_Scal
      real(real_kind),allocatable,   dimension(:,:)   :: All_centroid
      logical(log_kind),allocatable, dimension(:)     :: Hit_Good 
      integer(int_kind),allocatable, dimension(:,:)   :: All_Ngbrs_All, All_Ngbr
      real(real_kind),allocatable, dimension(:,:,:) :: Tot_Face_Centroid 

       ALLOCATE(All_Ngbrs_All(ncells_tot,nfc),         &
                All_Ngbr(ncells,nfc),                   &
                Tot_Vol(ncells_tot),                    &
                Tot_Scal(ncells_tot),                   &
                Hit_Good(ncells_tot),                   &
                All_centroid(ncells_tot,ndim),          &
                Tot_Face_Centroid(ncells_tot,ndim,nfc), & 
                STAT = memstat)
       call TLS_fatal_if_any ((memstat /= 0), 'KERN_CONVOLUTION_CENTER: Memory allocation error')

! ... Initialize arrays
       All_Ngbr = 0

! ... We need Ngbr_Cells_All of ALL cells in the domain
! ... fist pack on processor Ngbr_Cells_All into All_Ngbr, then collate & dist
! ... We'll build the ragged array later, right now punt

     do n=1,ncells
      do f=1,nfc
       All_Ngbr(n,f) = Mesh(n)%ngbr_Cell_orig(f)
      end do
     end do


    do f=1,nfc
     call PGSLIB_Collate(All_Ngbrs_All(:,f),All_Ngbr(:,f))
    end do
    do f=1,nfc
     call PGSLIB_Bcast(All_Ngbrs_All(:,f))
    enddo

! Now to Collate and Broadcast All Cell Centroids,Cell Volumes,Scalar Field

    call PGSLIB_Collate(Tot_Vol(:), Cell(:)%Volume)
    call PGSLIB_Bcast(Tot_Vol)

    call PGSLIB_Collate(Tot_Scal(:), Scal_Center(:))
    call PGSLIB_Bcast(Tot_Scal)

    do m = 1,ndim
      call PGSLIB_Collate(All_Centroid(:,m),Cell(:)%Centroid(m))
    end do
    do m = 1,ndim
      call PGSLIB_Bcast(All_Centroid(:,m))
    end do

    do f=1,nfc    
      do m=1,ndim
        call PGSLIB_Collate(Tot_Face_Centroid(:,m,f),Cell(:)%Face_Centroid(m,f))
      enddo
    enddo
    do f=1,nfc
      do m=1,ndim
        call PGSLIB_Bcast(Tot_Face_Centroid(:,m,f))
      enddo
    enddo  
 
! SWEEP finds the value of the convolution in each cell by using a
! recursive subroutine which examines cell neighbors and then the
! neighbors of neighbors

    SWEEP_VOF: do n=1,ncells

! Initialize Hit array and convolution value
! and set the centroid for the convolution
     Hit_Good = .false.
     KF_N = zero 
     KF_D = zero

     do m = 1,ndim
       cent(m) =  Cell(n)%Centroid(m)
     enddo 

! .. begin recursive loop

   RECURS_VOF:  do f = 1,nfc
          icell = Mesh(n)%Ngbr_Cell_orig(f)

          if(icell == DEGENERATE_FACE) cycle RECURS_VOF

            IF (icell == 0) THEN
!BC
              dist_x=2.*(Cell(n)%Face_Centroid(1,f)-cent(1))
              dist_y=2.*(Cell(n)%Face_Centroid(2,f)-cent(2))
              dist_z=2.*(Cell(n)%Face_Centroid(3,f)-cent(3))
              Ke   = Kern(dist_x,dist_y, dist_z,d)   
              KF_N = KF_N + Scal_Center(n)*Ke*Cell(n)%Volume
              KF_D = KF_D + Ke*Cell(n)%Volume

            ELSE
         
          if (icell > 0 .and. .not. Hit_Good(icell))then
            Hit_Good(icell) = .true.
            dist_x = all_centroid(icell,1) - cent(1)
            dist_y = all_centroid(icell,2) - cent(2)
            dist_z = all_centroid(icell,3) - cent(3)

            if (sqrt(dist_x**2+dist_y**2+dist_z**2) < d ) then             
              Ke   = Kern(dist_x,dist_y, dist_z,d)   
              KF_N = KF_N + Tot_Scal(icell)*Ke*Tot_vol(icell)
              KF_D = KF_D + Ke*Tot_vol(icell)
              call RECURS_CONV_CENTER(icell,Tot_Vol,Tot_Scal,All_Ngbrs_All,      &
                            All_centroid,Tot_Face_Centroid,cent,Hit_Good,KF_N,KF_D,d)
            endif
          endif

            END IF

          ENDDO RECURS_VOF
        
     if (KF_D.gt.zero) then
       Scal_Convol(n) = KF_N/KF_D
     else
       Scal_Convol(n)=zero
     endif 
! end loop on ncells
    End Do SWEEP_VOF


    DEALLOCATE(All_Ngbrs_All,      &
               All_Ngbr,           &
               Tot_Vol,            &
               Tot_Scal,           &
               Hit_Good,           &
               All_centroid,       &
               Tot_Face_Centroid)

  RETURN
  END SUBROUTINE KERN_CONVOLUTION_CENTER 

!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>
!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>

  RECURSIVE SUBROUTINE RECURS_CONV_CENTER(Cell_no,Tot_vol, Tot_scal,all_ngbr,    &
                                   all_cent,all_face_cent,cent,Hit_Array,KF_N,KF_D,d)
    !=======================================================================
    !  Purpose:
    !    Recursive subroutine to calculate the value of the convolution of 
    !    Tot_Scal with a kernel centered at some given cell -
    !=======================================================================
      use kind_module,       only: int_kind, real_kind, log_kind
      use mesh_module,       only: DEGENERATE_FACE
      Implicit none

      integer(int_kind),                  INTENT(IN) :: Cell_no
      real(real_kind),  dimension(:),     INTENT(IN) :: Tot_vol
      real(real_kind),  dimension(:),     INTENT(IN) :: Tot_Scal
      integer(int_kind),  dimension(:,:), INTENT(IN) :: all_ngbr
      real(real_kind),  dimension(:,:),   INTENT(IN) :: all_cent
      real(real_kind), dimension(:,:,:),  INTENT(IN) :: all_face_cent
      real(real_kind),  dimension(:),     INTENT(IN) :: cent
      logical(log_kind), dimension(:), INTENT(INOUT) :: Hit_Array
      real(real_kind),                 INTENT(INOUT) :: KF_N,KF_D
      real(real_kind),                    INTENT(IN) :: d 
 
! ... Local variables
      integer(int_kind)  :: icell,J
      real(real_kind)    :: radius,dist_x,dist_y,dist_z,Ke

! ... 
! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are


      RECURS_CONV_LOOP: Do J = 1,SIZE(all_ngbr(Cell_no,:))
         icell = all_ngbr(Cell_no,J)
         if(icell == DEGENERATE_FACE) cycle RECURS_CONV_LOOP

         if (icell==0) then
           dist_x = 2.*(all_face_cent(Cell_no,1,J) - all_cent(Cell_no,1))
           dist_y = 2.*(all_face_cent(Cell_no,2,J) - all_cent(Cell_no,2))
           dist_z = 2.*(all_face_cent(Cell_no,3,J) - all_cent(Cell_no,3))          
           Ke   = Kern(dist_x,dist_y, dist_z, d)
           KF_N = KF_N + Tot_Scal(Cell_no)*Ke*Tot_vol(Cell_no)
           KF_D = KF_D + Ke*Tot_vol(Cell_no)

         else

           if(icell > 0 .AND. .not. Hit_Array(icell))then
             Hit_Array(icell) = .true.
             dist_x = all_cent(icell,1) - cent(1)
             dist_y = all_cent(icell,2) - cent(2)
             dist_z = all_cent(icell,3) - cent(3)
             radius = sqrt(dist_x**2+dist_y**2+dist_z**2)
             if (radius < 1.05*d) then
               Ke   = Kern(dist_x,dist_y, dist_z, d)
               KF_N = KF_N + Tot_Scal(icell)*Ke*Tot_vol(icell)
               KF_D = KF_D + Ke*Tot_vol(icell)
               call RECURS_CONV_CENTER(icell,Tot_Vol,Tot_Scal,all_ngbr,all_cent,all_face_cent,cent,       &
                               Hit_Array,KF_N,KF_D,d)
             endif 
           endif 

         endif
 
        End Do RECURS_CONV_LOOP

      Return
  END SUBROUTINE RECURS_CONV_CENTER

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE KERN_INTERPOLATION_FACE (Scal_Center,Scal_Face,Scal_Flag,Face_Flag,d)

    !=======================================================================
    !  Purpose:
    !    Returns the values Phi at faces 
    !    interpolation from cell centers values using a kernel (Rudman 1998) 
    !=======================================================================
       use constants_module,        only: zero 
       use kind_module,             only: real_kind, int_kind, log_kind
       use mesh_module,             only: Mesh,Cell, degenerate_face
       use pgslib_module,           only: PGSLIB_BCast, PGSLIB_Collate
       use parameter_module,        only: ncells, ncells_tot, ndim, nfc
       use var_vector_module

      IMPLICIT NONE

       real(real_kind), dimension(ncells),          INTENT(IN)    :: Scal_Center
       real(real_kind), dimension(nfc,ncells),      INTENT(INOUT) :: Scal_Face 
       logical(log_kind), dimension(ncells),        INTENT(IN)    :: Scal_Flag  
       logical(log_kind),dimension(nfc,ncells),     INTENT(IN)    :: Face_Flag
       real(real_kind)                                            :: d

! local variables
      integer(int_kind)                    :: n,f,fo,m,icell, memstat
      real(real_kind)                      :: dist_x, dist_y, dist_z, Ke,KF_N, KF_D 
      real(real_kind),dimension(ndim,nfc)  :: face_cent
      real(real_kind),dimension(ndim)      :: fcent

      real(real_kind),allocatable,   dimension(:)   :: Tot_Vol, Tot_Scal
      real(real_kind),allocatable,   dimension(:,:) :: All_centroid
      logical(log_kind),allocatable, dimension(:)   :: Hit_Good 
      integer(int_kind),allocatable, dimension(:,:) :: All_Ngbrs_All, All_Ngbr
      logical(log_kind),allocatable, dimension(:)    :: All_Scal_Flag



       ALLOCATE(All_Ngbrs_All(ncells_tot,nfc),   &
                All_Ngbr(ncells,nfc),            &
                Tot_Vol(ncells_tot),             &
                Tot_Scal(ncells_tot),            &
                Hit_Good(ncells_tot),            &
                All_centroid(ncells_tot,ndim),   &
                All_Scal_Flag(ncells_tot),       & 
                STAT = memstat)
       call TLS_fatal_if_any ((memstat /= 0), 'KERN_INTERPOLATION_FACE: Memory allocation error')

! ... Initialize arrays
       All_Ngbr = 0

! ... We need Ngbr_Cells_All of ALL cells in the domain
! ... fist pack on processor Ngbr_Cells_All into All_Ngbr, then collate & dist
! ... We'll build the ragged array later, right now punt

     do n=1,ncells
      do f=1,nfc
       All_Ngbr(n,f) = Mesh(n)%ngbr_Cell_orig(f)
      end do
     end do


    do f=1,nfc
     call PGSLIB_Collate(All_Ngbrs_All(:,f),All_Ngbr(:,f))
    end do
    do f=1,nfc
     call PGSLIB_Bcast(All_Ngbrs_All(:,f))
    enddo

! Now to Collate and Broadcast All Cell Centroids,Cell Volumes,Scalar Field

    call PGSLIB_Collate(Tot_Vol(:), Cell(:)%Volume)
    call PGSLIB_Bcast(Tot_Vol)

    call PGSLIB_Collate(Tot_Scal(:), Scal_Center(:))
    call PGSLIB_Bcast(Tot_Scal)

    do m = 1,ndim
      call PGSLIB_Collate(All_Centroid(:,m),Cell(:)%Centroid(m))
    end do
    do m = 1,ndim
      call PGSLIB_Bcast(All_Centroid(:,m))
    end do
      
    call PGSLIB_Collate(All_Scal_Flag(:),Scal_Flag(:))
    call PGSLIB_Bcast(All_Scal_Flag(:))
 
! SWEEP finds the value of the convolution in each cell by using a
! recursive subroutine which examines cell neighbors and then the
! neighbors of neighbors

    SWEEP_VOF: do n=1,ncells

! outer loop on faces
     FACE_LOOP: do fo=1,nfc
        Scal_Face(fo,n)=zero
        if (Face_Flag(fo,n)) then
          Hit_Good = .false.
          KF_N = zero
          KF_D = zero
          do m=1,ndim
            face_cent(m,fo)=Cell(n)%Face_Centroid(m,fo)
          enddo

!first consider current cell
          Hit_Good(n)=.true.
          fcent(:)=face_cent(:,fo)
          dist_x = All_centroid(n,1) - fcent(1)
          dist_y = All_centroid(n,2) - fcent(2)
          dist_z = All_centroid(n,3) - fcent(3)           

          if (sqrt(dist_x**2+dist_y**2+dist_z**2) < d ) then 
            if (All_Scal_Flag(n)) then
              Ke   = Kern(dist_x,dist_y, dist_z,d)
              KF_N = KF_N + Tot_Scal(n)*Ke*Tot_vol(n)
              KF_D = KF_D + Ke*Tot_vol(n)
            endif
          endif
 
! .. begin recursive loop
   RECURS_VOF:  do f = 1,nfc
          icell = Mesh(n)%Ngbr_Cell_orig(f)
          if(icell == DEGENERATE_FACE) cycle RECURS_VOF
         
          if (icell > 0 .and. .not. Hit_Good(icell))then
            Hit_Good(icell) = .true.
            fcent(:)=face_cent(:,fo)
            dist_x = All_centroid(icell,1) - fcent(1)
            dist_y = All_centroid(icell,2) - fcent(2)
            dist_z = All_centroid(icell,3) - fcent(3)

            if (sqrt(dist_x**2+dist_y**2+dist_z**2) < d) then
              if (All_Scal_Flag(icell)) then
                Ke   = Kern(dist_x,dist_y, dist_z,d)   
                KF_N = KF_N + Tot_Scal(icell)*Ke*Tot_vol(icell)
                KF_D = KF_D + Ke*Tot_vol(icell) 
              endif 
              call RECURS_CONV_FACE(icell,Tot_Vol,Tot_Scal,All_Ngbrs_All,      &
                               All_centroid,fcent,Hit_Good,KF_N,KF_D,d,        &          
                               All_Scal_Flag)
            endif
          endif

!            END IF

          ENDDO RECURS_VOF
          if (KF_D.gt.zero) then 
            Scal_Face(fo,n) = KF_N/KF_D
          else
            Scal_Face(fo,n)=zero
          endif
        endif 
! end outer loop on faces
      end do FACE_LOOP

! end loop on ncells
    End Do SWEEP_VOF


    DEALLOCATE(All_Ngbrs_All,      &
               All_Ngbr,           &
               Tot_Vol,            &
               Tot_Scal,           &
               Hit_Good,           &
               All_centroid,       &
               ALL_Scal_Flag)

  RETURN
  END SUBROUTINE KERN_INTERPOLATION_FACE 

!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>

  RECURSIVE SUBROUTINE RECURS_CONV_FACE(Cell_no,Tot_vol, Tot_scal,all_ngbr,    &
                                   all_cent,fcent,Hit_Array,KF_N,KF_D,d,     &
                                   All_Scal_Flag)
    !=======================================================================
    !  Purpose:
    !    Recursive subroutine to calculate the value of the convolution of 
    !    Tot_Scal with a kernel centered at some given cell -
    !=======================================================================
    use kind_module,       only: int_kind, real_kind, log_kind
    use mesh_module,       only: DEGENERATE_FACE
    Implicit none

      integer(int_kind),                  INTENT(IN) :: Cell_no
      real(real_kind),  dimension(:),     INTENT(IN) :: Tot_vol
      real(real_kind),  dimension(:),     INTENT(IN) :: Tot_Scal
      integer(int_kind),  dimension(:,:), INTENT(IN) :: all_ngbr
      real(real_kind),  dimension(:,:),   INTENT(IN) :: all_cent
      real(real_kind),  dimension(:),     INTENT(IN) :: fcent
      logical(log_kind), dimension(:), INTENT(INOUT) :: Hit_Array
      real(real_kind),                 INTENT(INOUT) :: KF_N,KF_D
      real(real_kind),                    INTENT(IN) :: d 
      logical(log_kind), dimension(:),    INTENT(IN) :: All_Scal_Flag
 
! ... Local variables
      integer(int_kind)  :: icell,J
      real(real_kind)    :: radius,dist_x,dist_y,dist_z,Ke

! ... 
! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are


      RECURS_CONV_LOOP: Do J = 1,SIZE(all_ngbr(Cell_no,:))
         icell = all_ngbr(Cell_no,J)
         if(icell == DEGENERATE_FACE) cycle RECURS_CONV_LOOP

          if(icell > 0 .AND. .not. Hit_Array(icell))then
            Hit_Array(icell) = .true.
            dist_x = all_cent(icell,1) - fcent(1)
            dist_y = all_cent(icell,2) - fcent(2)
            dist_z = all_cent(icell,3) - fcent(3)
            radius = sqrt(dist_x**2+dist_y**2+dist_z**2)
            if (radius < 1.05*d) then
              if (All_Scal_Flag(icell)) then
                Ke   = Kern(dist_x,dist_y, dist_z, d)
                KF_N = KF_N + Tot_Scal(icell)*Ke*Tot_vol(icell)
                KF_D = KF_D + Ke*Tot_vol(icell)
              endif  
              call RECURS_CONV_FACE(icell,Tot_Vol,Tot_Scal,all_ngbr,all_cent,fcent, &
                                    Hit_Array,KF_N,KF_D,d,All_Scal_Flag)
            endif 
          endif 

        End Do RECURS_CONV_LOOP

      Return
  END SUBROUTINE RECURS_CONV_FACE

!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>

      real FUNCTION Kern(x,y,z,d)
!Function that returns that computes the Kernel values
!kernel as given in Rudman IJNMF 1998

      use kind_module,       only: real_kind
      use constants_module,  only: pi,zero
      implicit none

      real(real_kind) :: x,y,z,d
      real(real_kind) :: rdist
      

      rdist=sqrt(x**2+y**2+z**2)
      if (rdist.lt.d/2.) then
        Kern=40./(7.*pi)*(1-6.*(rdist/d)**2+6.*(rdist/d)**3)
      endif
      if (rdist.ge.d/2.and.rdist.lt.d) then
        Kern=80./(7.*pi)*(1-(rdist/d))**3
      endif
      if (rdist.ge.d) then
        Kern=zero
      endif
      Kern=Kern/(d**2)

      END FUNCTION Kern


!<><><><><><>><<><>><><><><>><><><><><><><><><><><><><><><><>><><><><><<><>

END MODULE KERNEL_INTERPOLATION_MODULE

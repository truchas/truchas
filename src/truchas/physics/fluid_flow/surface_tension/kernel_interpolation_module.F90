!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kernel_interpolation_module 
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
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: kern_convolution_center, kern_interpolation_face 

contains

  !=======================================================================
  !  Purpose:
  !    Returns the convolved value of Scal_Center
  !    using Rudman 1998 kernel 
  !=======================================================================

  subroutine kern_convolution_center (Scal_Center, Scal_Convol, d)

    use mesh_module, only: Mesh,Cell, degenerate_face
    use pgslib_module, only: PGSLIB_BCast, PGSLIB_Collate
    use parameter_module, only: ncells, ncells_tot, ndim, nfc
    use var_vector_module

    real(r8), dimension(ncells), intent(in)    :: Scal_Center
    real(r8), dimension(ncells), intent(inout) :: Scal_Convol 
    real(r8) :: d
    
    integer  :: n, f, m, icell, memstat
    real(r8) :: dist_x, dist_y, dist_z, Ke, KF_N, KF_D, cent(ndim)
    real(r8), allocatable :: Tot_Vol(:), Tot_Scal(:), All_centroid(:,:), Tot_Face_Centroid(:,:,:)
    logical,  allocatable :: Hit_Good(:)
    integer,  allocatable :: All_Ngbrs_All(:,:), All_Ngbr(:,:)

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

    do n=1, ncells
      do f=1, nfc
        All_Ngbr(n,f) = Mesh(n)%ngbr_Cell_orig(f)
      end do
    end do


    do f=1,nfc
      call PGSLIB_Collate(All_Ngbrs_All(:,f),All_Ngbr(:,f))
    end do
    do f=1,nfc
      call PGSLIB_Bcast(All_Ngbrs_All(:,f))
    end do

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
      end do
    end do
    do f=1,nfc
      do m=1,ndim
        call PGSLIB_Bcast(Tot_Face_Centroid(:,m,f))
      end do
    end do  
 
! SWEEP finds the value of the convolution in each cell by using a
! recursive subroutine which examines cell neighbors and then the
! neighbors of neighbors

    SWEEP_VOF: do n=1, ncells

! Initialize Hit array and convolution value
! and set the centroid for the convolution
      Hit_Good = .false.
      KF_N = 0.0_r8 
      KF_D = 0.0_r8

      do m = 1,ndim
        cent(m) =  Cell(n)%Centroid(m)
      end do 

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

          if (icell > 0 .and. .not.Hit_Good(icell)) then
            Hit_Good(icell) = .true.
            dist_x = all_centroid(icell,1) - cent(1)
            dist_y = all_centroid(icell,2) - cent(2)
            dist_z = all_centroid(icell,3) - cent(3)
            if (sqrt(dist_x**2+dist_y**2+dist_z**2) < d) then             
              Ke   = Kern(dist_x,dist_y, dist_z,d)   
              KF_N = KF_N + Tot_Scal(icell)*Ke*Tot_vol(icell)
              KF_D = KF_D + Ke*Tot_vol(icell)
              call RECURS_CONV_CENTER(icell,Tot_Vol,Tot_Scal,All_Ngbrs_All,      &
                            All_centroid,Tot_Face_Centroid,cent,Hit_Good,KF_N,KF_D,d)
            endif
          endif
        end if

      end do RECURS_VOF

     if (KF_D > 0.0_r8) then
       Scal_Convol(n) = KF_N/KF_D
     else
       Scal_Convol(n)=0.0_r8
     endif 
     
    end do SWEEP_VOF

    deallocate(All_Ngbrs_All,      &
               All_Ngbr,           &
               Tot_Vol,            &
               Tot_Scal,           &
               Hit_Good,           &
               All_centroid,       &
               Tot_Face_Centroid)

  end subroutine kern_convolution_center 

  !=======================================================================
  !  Purpose:
  !    Recursive subroutine to calculate the value of the convolution of 
  !    Tot_Scal with a kernel centered at some given cell -
  !=======================================================================

  recursive subroutine recurs_conv_center (Cell_no,Tot_vol, Tot_scal,all_ngbr,    &
                                   all_cent,all_face_cent,cent,Hit_Array,KF_N,KF_D,d)
    use mesh_module, only: DEGENERATE_FACE

    integer,  intent(in) :: Cell_no
    real(r8), intent(in) :: Tot_vol(:)
    real(r8), intent(in) :: Tot_Scal(:)
    integer,  intent(in) :: all_ngbr(:,:)
    real(r8), intent(in) :: all_cent(:,:)
    real(r8), intent(in) :: all_face_cent(:,:,:)
    real(r8), intent(in) :: cent(:)
    logical,  intent(inout) :: Hit_Array(:)
    real(r8), intent(inout) :: KF_N, KF_D
    real(r8), intent(in) :: d 

! ... Local variables
    integer  :: icell,J
    real(r8) :: radius,dist_x,dist_y,dist_z,Ke

! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are

    RECURS_CONV_LOOP: do j = 1, size(all_ngbr(Cell_no,:))
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
          end if 
        end if 
      end if
    end do RECURS_CONV_LOOP

  end subroutine recurs_conv_center

  !=======================================================================
  !  Purpose:
  !    Returns the values Phi at faces 
  !    interpolation from cell centers values using a kernel (Rudman 1998) 
  !=======================================================================

  subroutine kern_interpolation_face (Scal_Center,Scal_Face,Scal_Flag,Face_Flag,d)

    use mesh_module, only: Mesh,Cell, degenerate_face
    use pgslib_module, only: PGSLIB_BCast, PGSLIB_Collate
    use parameter_module, only: ncells, ncells_tot, ndim, nfc
    use var_vector_module

    real(r8), dimension(ncells), INTENT(IN) :: Scal_Center
    real(r8), dimension(nfc,ncells), INTENT(INOUT) :: Scal_Face 
    logical, dimension(ncells), INTENT(IN) :: Scal_Flag  
    logical, dimension(nfc,ncells), INTENT(IN) :: Face_Flag
    real(r8) :: d

! local variables
    integer :: n,f,fo,m,icell, memstat
    real(r8) :: dist_x, dist_y, dist_z, Ke,KF_N, KF_D 
    real(r8), dimension(ndim,nfc) :: face_cent
    real(r8), dimension(ndim) :: fcent

    real(r8), allocatable :: Tot_Vol(:), Tot_Scal(:), All_centroid(:,:)
    logical,  allocatable :: Hit_Good(:), All_Scal_Flag(:)
    integer,  allocatable :: All_Ngbrs_All(:,:), All_Ngbr(:,:)

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
    end do

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
        Scal_Face(fo,n)=0.0_r8
        if (Face_Flag(fo,n)) then
          Hit_Good = .false.
          KF_N = 0.0_r8
          KF_D = 0.0_r8
          do m=1,ndim
            face_cent(m,fo)=Cell(n)%Face_Centroid(m,fo)
          end do

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
            end if
          end if
 
! .. begin recursive loop
          RECURS_VOF: do f=1, nfc
            icell = Mesh(n)%Ngbr_Cell_orig(f)
            if (icell == DEGENERATE_FACE) cycle RECURS_VOF

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
                end if 
                call RECURS_CONV_FACE(icell,Tot_Vol,Tot_Scal,All_Ngbrs_All,      &
                                 All_centroid,fcent,Hit_Good,KF_N,KF_D,d,        &          
                                 All_Scal_Flag)
              end if
            end if
          end do RECURS_VOF
          if (KF_D > 0.0_r8) then 
            Scal_Face(fo,n) = KF_N/KF_D
          else
            Scal_Face(fo,n)=0.0_r8
          end if
        end if 
      end do FACE_LOOP
    end do SWEEP_VOF

    DEALLOCATE(All_Ngbrs_All,      &
               All_Ngbr,           &
               Tot_Vol,            &
               Tot_Scal,           &
               Hit_Good,           &
               All_centroid,       &
               ALL_Scal_Flag)

  end subroutine kern_interpolation_face 

  !=======================================================================
  !  Purpose:
  !    Recursive subroutine to calculate the value of the convolution of 
  !    Tot_Scal with a kernel centered at some given cell -
  !=======================================================================

  recursive subroutine recurs_conv_face (Cell_no,Tot_vol, Tot_scal,all_ngbr,    &
                                   all_cent,fcent,Hit_Array,KF_N,KF_D,d,     &
                                   All_Scal_Flag)
    use mesh_module, only: DEGENERATE_FACE

    integer,  intent(in) :: Cell_no
    real(r8), intent(in) :: Tot_vol(:)
    real(r8), intent(in) :: Tot_Scal(:)
    integer,  intent(in) :: all_ngbr(:,:)
    real(r8), intent(in) :: all_cent(:,:)
    real(r8), intent(in) :: fcent(:)
    logical,  intent(inout) :: Hit_Array(:)
    real(r8), intent(inout) :: KF_N,KF_D
    real(r8), intent(in) :: d 
    logical,  intent(in) :: All_Scal_Flag(:)

! ... Local variables
    integer  :: icell,J
    real(r8) :: radius, dist_x, dist_y, dist_z, Ke

! ... check neighbor cells of Cell_no, if they are in the domain and 
! ... haven't been hit check to see if they are within the support
! ... continue recursion if they are

    RECURS_CONV_LOOP: do j = 1, size(all_ngbr(Cell_no,:))
      icell = all_ngbr(Cell_no,J)
      if(icell == DEGENERATE_FACE) cycle RECURS_CONV_LOOP
      if(icell > 0 .and. .not.Hit_Array(icell)) then
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
          end if  
          call RECURS_CONV_FACE(icell,Tot_Vol,Tot_Scal,all_ngbr,all_cent,fcent, &
                                Hit_Array,KF_N,KF_D,d,All_Scal_Flag)
        end if 
      end if 
    end do RECURS_CONV_LOOP

  end subroutine recurs_conv_face

  ! Function that returns that computes the Kernel values
  ! kernel as given in Rudman IJNMF 1998

  real function Kern (x, y, z, d)

    use constants_module, only: pi

    real(r8) :: x ,y, z, d
    real(r8) :: rdist

    rdist=sqrt(x**2+y**2+z**2)
    if (rdist.lt.d/2.) then
      Kern=40./(7.*pi)*(1-6.*(rdist/d)**2+6.*(rdist/d)**3)
    end if
    if (rdist.ge.d/2.and.rdist.lt.d) then
      Kern=80./(7.*pi)*(1-(rdist/d))**3
    end if
    if (rdist.ge.d) then
      Kern=0.0_r8
    end if
    Kern=Kern/(d**2)

  end function Kern

end module kernel_interpolation_module

!!
!! INDEX_MAP_TYPE
!!
!! This module provides core capabilities that support the use of distributed
!! arrays in MPI-based SPMD programs through the INDEX_MAP derived type that
!! describes the mapping of an array's index set to processes. The mapping
!! allows for overlap between processes, and provides collective gather and
!! scatter procedures associated with that overlap. The module provides
!! additional procedures for distributing and localizing serial indirect
!! indexing arrays.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

module index_map_type

  use,intrinsic :: iso_fortran_env, only: i4 => int32, i8 => int64, r4 => real32, r8 => real64
  use mpi
  implicit none
  private

  type, public :: index_map
    integer :: comm = MPI_COMM_WORLD !TODO: passed to INIT
    integer :: nproc, rank, root=0
    logical, private :: is_root
    integer :: onp_size=0, offp_size=0, local_size=0  ! PROTECTED
    integer :: global_size=0, first_gid=0, last_gid=0 ! PROTECTED; int64 option?
    ! off-process gather/scatter communication data
    integer, private :: gather_comm, scatter_comm
    integer, allocatable :: offp_index(:)  ! int64 option? !TODO: private?
    integer, allocatable, private :: offp_counts(:), offp_displs(:)
    integer, allocatable, private :: onp_index(:), onp_counts(:), onp_displs(:)
    ! gather/scatter communication data
    integer, allocatable, private :: counts(:), displs(:)
  contains
    generic   :: init => init_dist, init_root, init_dist_offp, init_root_offp, init_ragged
    procedure :: add_offp_index
    procedure :: global_index
    procedure :: defined
    generic :: gather_offp => &
        gath1_i4_1, gath2_i4_1, gath1_i4_2, gath2_i4_2, gath1_i4_3, gath2_i4_3, &
        gath1_r4_1, gath2_r4_1, gath1_r4_2, gath2_r4_2, gath1_r4_3, gath2_r4_3, &
        gath1_r8_1, gath2_r8_1, gath1_r8_2, gath2_r8_2, gath1_r8_3, gath2_r8_3, &
        gath1_c4_1, gath2_c4_1, gath1_c4_2, gath2_c4_2, gath1_c4_3, gath2_c4_3, &
        gath1_c8_1, gath2_c8_1, gath1_c8_2, gath2_c8_2, gath1_c8_3, gath2_c8_3, &
        gath1_dl_1, gath2_dl_1, gath1_dl_2, gath2_dl_2, gath1_dl_3, gath2_dl_3
    generic :: scatter_offp_sum => &
        scat1_sum_i4_1, scat2_sum_i4_1, &
        scat1_sum_r4_1, scat2_sum_r4_1, &
        scat1_sum_r8_1, scat2_sum_r8_1, &
        scat1_sum_c4_1, scat2_sum_c4_1, &
        scat1_sum_c8_1, scat2_sum_c8_1
    generic :: scatter_offp_min => &
        scat1_min_i4_1, scat2_min_i4_1, &
        scat1_min_r4_1, scat2_min_r4_1, &
        scat1_min_r8_1, scat2_min_r8_1
    generic :: scatter_offp_max => &
        scat1_max_i4_1, scat2_max_i4_1, &
        scat1_max_r4_1, scat2_max_r4_1, &
        scat1_max_r8_1, scat2_max_r8_1
    generic :: scatter_offp_or => scat1_or_dl_1, scat2_or_dl_1
    generic :: scatter_offp_and => scat1_and_dl_1, scat2_and_dl_1
    generic :: scatter => &
        scat_i4_1, scat_i8_1, scat_r4_1, scat_r8_1, scat_c4_1, scat_c8_1, scat_dl_1, &
        scat_i4_2, scat_i8_2, scat_r4_2, scat_r8_2, scat_c4_2, scat_c8_2, scat_dl_2, &
        scat_i4_3, scat_i8_3, scat_r4_3, scat_r8_3, scat_c4_3, scat_c8_3, scat_dl_3
    generic :: gather => &
        gath_i4_1, gath_i8_1, gath_r4_1, gath_r8_1, gath_c4_1, gath_c8_1, gath_dl_1, &
        gath_i4_2, gath_i8_2, gath_r4_2, gath_r8_2, gath_c4_2, gath_c8_2, gath_dl_2, &
        gath_i4_3, gath_i8_3, gath_r4_3, gath_r8_3, gath_c4_3, gath_c8_3, gath_dl_3
    generic :: localize_index_array => localize_index_array_serial_1, localize_index_array_serial_2, &
        localize_index_array_dist_1, localize_index_array_dist_2, localize_index_struct_serial
    procedure, private :: init_dist, init_root, init_dist_offp, init_root_offp, init_ragged
    procedure, private :: &
        gath1_i4_1, gath2_i4_1, gath1_i4_2, gath2_i4_2, gath1_i4_3, gath2_i4_3, &
        gath1_r4_1, gath2_r4_1, gath1_r4_2, gath2_r4_2, gath1_r4_3, gath2_r4_3, &
        gath1_r8_1, gath2_r8_1, gath1_r8_2, gath2_r8_2, gath1_r8_3, gath2_r8_3, &
        gath1_c4_1, gath2_c4_1, gath1_c4_2, gath2_c4_2, gath1_c4_3, gath2_c4_3, &
        gath1_c8_1, gath2_c8_1, gath1_c8_2, gath2_c8_2, gath1_c8_3, gath2_c8_3, &
        gath1_dl_1, gath2_dl_1, gath1_dl_2, gath2_dl_2, gath1_dl_3, gath2_dl_3
    procedure, private :: &
        scat1_sum_i4_1, scat2_sum_i4_1, &
        scat1_sum_r4_1, scat2_sum_r4_1, &
        scat1_sum_r8_1, scat2_sum_r8_1, &
        scat1_sum_c4_1, scat2_sum_c4_1, &
        scat1_sum_c8_1, scat2_sum_c8_1
    procedure, private :: &
        scat1_min_i4_1, scat2_min_i4_1, &
        scat1_min_r4_1, scat2_min_r4_1, &
        scat1_min_r8_1, scat2_min_r8_1
    procedure, private :: &
        scat1_max_i4_1, scat2_max_i4_1, &
        scat1_max_r4_1, scat2_max_r4_1, &
        scat1_max_r8_1, scat2_max_r8_1
    procedure, private :: scat1_or_dl_1, scat2_or_dl_1
    procedure, private :: scat1_and_dl_1, scat2_and_dl_1
    procedure, private :: &
        scat_i4_1, scat_i8_1, scat_r4_1, scat_r8_1, scat_c4_1, scat_c8_1, scat_dl_1, &
        scat_i4_2, scat_i8_2, scat_r4_2, scat_r8_2, scat_c4_2, scat_c8_2, scat_dl_2, &
        scat_i4_3, scat_i8_3, scat_r4_3, scat_r8_3, scat_c4_3, scat_c8_3, scat_dl_3
    procedure, private :: &
        gath_i4_1, gath_i8_1, gath_r4_1, gath_r8_1, gath_c4_1, gath_c8_1, gath_dl_1, &
        gath_i4_2, gath_i8_2, gath_r4_2, gath_r8_2, gath_c4_2, gath_c8_2, gath_dl_2, &
        gath_i4_3, gath_i8_3, gath_r4_3, gath_r8_3, gath_c4_3, gath_c8_3, gath_dl_3
    procedure, private :: localize_index_array_serial_1, localize_index_array_serial_2, &
        localize_index_array_dist_1, localize_index_array_dist_2, localize_index_struct_serial
    procedure, private :: add_offp_index_set ! Type bound to workaound gfortran bug
  end type

  interface
    module subroutine gath1_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_c4_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_c8_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine

    module subroutine gath2_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:)
      integer(i4), intent(inout) :: offp_data(:)
    end subroutine
    module subroutine gath2_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:)
      real(r4), intent(inout) :: offp_data(:)
    end subroutine
    module subroutine gath2_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:)
      real(r8), intent(inout) :: offp_data(:)
    end subroutine
    module subroutine gath2_c4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:)
      complex(r4), intent(inout) :: offp_data(:)
    end subroutine
    module subroutine gath2_c8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:)
      complex(r8), intent(inout) :: offp_data(:)
    end subroutine
    module subroutine gath2_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:)
      logical, intent(inout) :: offp_data(:)
    end subroutine

    module subroutine gath1_i4_2(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_r4_2(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_r8_2(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_c4_2(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_c8_2(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_dl_2(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:,:)
    end subroutine

    module subroutine gath2_i4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:,:)
      integer(i4), intent(inout) :: offp_data(:,:)
    end subroutine
    module subroutine gath2_r4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:,:)
      real(r4), intent(inout) :: offp_data(:,:)
    end subroutine
    module subroutine gath2_r8_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:,:)
      real(r8), intent(inout) :: offp_data(:,:)
    end subroutine
    module subroutine gath2_c4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:,:)
      complex(r4), intent(inout) :: offp_data(:,:)
    end subroutine
    module subroutine gath2_c8_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:,:)
      complex(r8), intent(inout) :: offp_data(:,:)
    end subroutine
    module subroutine gath2_dl_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:,:)
      logical, intent(inout) :: offp_data(:,:)
    end subroutine

    module subroutine gath1_i4_3(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_r4_3(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_r8_3(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_c4_3(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_c8_3(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_dl_3(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:,:,:)
    end subroutine

    module subroutine gath2_i4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:,:,:)
      integer(i4), intent(inout) :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_r4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:,:,:)
      real(r4), intent(inout) :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_r8_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:,:,:)
      real(r8), intent(inout) :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_c4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:,:,:)
      complex(r4), intent(inout) :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_c8_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:,:,:)
      complex(r8), intent(inout) :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_dl_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:,:,:)
      logical, intent(inout) :: offp_data(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine scat1_sum_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_c4_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_c4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: onp_data(:)
      complex(r4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_c8_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_c8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: onp_data(:)
      complex(r8), intent(in) :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_min_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_min_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_min_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in) :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_max_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_max_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scat1_max_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in) :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_or_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_or_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in) :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_and_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_and_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in) :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat_i4_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:)
      integer(i4), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_i8_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:)
      integer(i8), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_r4_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:)
      real(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_r8_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:)
      real(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_c4_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:)
      complex(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_c8_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:)
      complex(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_dl_1(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:)
      logical, intent(inout) :: dest(:)
    end subroutine
    module subroutine scat_i4_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:)
      integer(i4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_i8_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:)
      integer(i8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_r4_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:)
      real(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_r8_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:)
      real(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_c4_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:)
      complex(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_c8_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:)
      complex(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_dl_2(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:)
      logical, intent(inout) :: dest(:,:)
    end subroutine
    module subroutine scat_i4_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:,:)
      integer(i4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_i8_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:,:)
      integer(i8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_r4_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:,:)
      real(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_r8_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:,:)
      real(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_c4_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:,:)
      complex(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_c8_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:,:)
      complex(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine scat_dl_3(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:,:)
      logical, intent(inout) :: dest(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine gath_i4_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:)
      integer(i4), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_i8_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:)
      integer(i8), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_r4_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:)
      real(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_r8_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:)
      real(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_c4_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:)
      complex(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_c8_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:)
      complex(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_dl_1(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:)
      logical, intent(inout) :: dest(:)
    end subroutine
    module subroutine gath_i4_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:)
      integer(i4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_i8_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:)
      integer(i8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_r4_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:)
      real(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_r8_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:)
      real(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_c4_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:)
      complex(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_c8_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:)
      complex(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_dl_2(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:)
      logical, intent(inout) :: dest(:,:)
    end subroutine
    module subroutine gath_i4_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:,:)
      integer(i4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_i8_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:,:)
      integer(i8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_r4_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:,:)
      real(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_r8_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:,:)
      real(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_c4_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:,:)
      complex(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_c8_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:,:)
      complex(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine gath_dl_3(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:,:)
      logical, intent(inout) :: dest(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine localize_index_array_serial_1(domain, g_index, range, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:)
      integer, allocatable, intent(out) :: l_index(:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_serial_2(domain, g_index, range, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:,:)
      integer, allocatable, intent(out) :: l_index(:,:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_dist_1(range, index, stat)
      class(index_map), intent(inout) :: range
      integer, intent(inout) :: index(:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_dist_2(range, index, stat)
      class(index_map), intent(inout) :: range
      integer, contiguous, intent(inout), target :: index(:,:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_struct_serial(domain, g_count, g_index, range, l_count, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:), g_count(:)
      integer, allocatable, intent(out) :: l_index(:), l_count(:)
      integer, intent(out), optional :: stat
    end subroutine
  end interface

contains

  !! Each rank supplied its block size
  subroutine init_dist(this, bsize, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize
    integer, intent(in), optional :: root

    integer :: ierr

    !this%comm = comm -- default initialized to MPI_COMM_WORLD for now
    call MPI_Comm_size(this%comm, this%nproc, ierr)
    call MPI_Comm_rank(this%comm, this%rank, ierr)
    if (present(root)) this%root = root ! overwrite default
    this%is_root = (this%rank == this%root)

    this%onp_size = bsize
    this%offp_size = 0
    this%local_size = this%onp_size + this%offp_size
    call MPI_Scan(bsize, this%last_gid, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
    this%first_gid = this%last_gid - this%onp_size + 1
    this%global_size = this%last_gid
    call MPI_Bcast(this%global_size, 1, MPI_INTEGER, this%nproc-1, this%comm, ierr)

    call add_dist_coll_info(this)

  end subroutine

  !! One root rank has an array of rank block sizes
  subroutine init_root(this, bsizes, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsizes(:)
    integer, intent(in), optional :: root
    integer :: ierr, rank, nproc, bsize
    if (present(root)) this%root = root ! overwrite default
    call MPI_Comm_rank(this%comm, rank, ierr)
    if (rank == this%root) then
      call MPI_Comm_size(this%comm, nproc, ierr)
      INSIST(size(bsizes) == nproc)
    end if
    call MPI_Scatter(bsizes, 1, MPI_INTEGER, bsize, 1, MPI_INTEGER, this%root, this%comm, ierr)
    call init_dist(this, bsize, root)
  end subroutine

  !! Each rank supplied with its block size and list of off-process indices
  subroutine init_dist_offp(this, bsize, offp_index, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize
    integer, intent(in) :: offp_index(:) ! 64-bit option?
    integer, intent(in), optional :: root
    call init_dist(this, bsize, root)
    call add_offp_index(this, offp_index)
  end subroutine

  !! One root rank has an array of rank block sizes and off-process indices for all ranks
  subroutine init_root_offp(this, bsizes, offp_counts, offp_indices, root)

    class(index_map), intent(out) :: this
    integer, intent(in) :: bsizes(:), offp_counts(:)
    integer, intent(in) :: offp_indices(:)
    integer, intent(in), optional :: root
    type(index_map) :: temp
    integer, allocatable :: offp_index(:)
    call init_root(this, bsizes, root)
    call init_root(temp, offp_counts, root)
    allocate(offp_index(temp%onp_size))
    call temp%scatter(offp_indices, offp_index)
    call add_offp_index(this, offp_index)
  end subroutine

  !! Ragged index set based on another index set
  subroutine init_ragged(this, domain, g_count)

    class(index_map), intent(out) :: this
    type(index_map), intent(in) :: domain
    integer, intent(in) :: g_count(:)

    integer :: n, i, j, ierr, nmax
    integer, allocatable :: l_count(:), offset(:), offp_index(:)

    ASSERT(size(g_count) >= merge(domain%global_size,0,domain%is_root))

    allocate(l_count(domain%local_size))
    call domain%scatter(g_count, l_count)
    if (allocated(domain%offp_index)) call domain%gather_offp(l_count)

    call this%init(sum(l_count(1:domain%onp_size)))

    if (allocated(domain%offp_index)) then
      call MPI_Scan(this%onp_size, n, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
      n = n - this%onp_size ! exclusive scan of this%onp_size
      allocate(offset(domain%local_size))
      if (size(offset) > 0) then
        offset(1) = n
        do j = 2, domain%onp_size
          offset(j) = offset(j-1) + l_count(j-1)
        end do
        call domain%gather_offp(offset)
      end if
      n = sum(l_count(domain%onp_size+1:domain%local_size))
      call MPI_Allreduce(n, nmax, 1, MPI_INTEGER, MPI_MAX, domain%comm, ierr)
      if (nmax > 0) then ! have off-process indices for this index map
        allocate(offp_index(n))
        n = 0
        do j = domain%onp_size + 1, domain%local_size
          offp_index(n+1:n+l_count(j)) = [(offset(j)+i, i=1,l_count(j))]
          n = n + l_count(j)
        end do
        call add_offp_index(this, offp_index)
      end if
    end if

  end subroutine

  !! Return the global index corresponding to local index N
  elemental function global_index(this, n) result(gid)
    class(index_map), intent(in) :: this
    integer, intent(in) :: n
    integer :: gid
    gid = -1
    if (n < 1) return
    if (n <= this%onp_size) then
      gid = this%first_gid + n - 1
    else if (n <= this%local_size) then
      gid = this%offp_index(n-this%onp_size)
    end if
  end function

  !! Add off-process indices to an already initialized index map.
  !! This is a public interface with untrusted OFFP_INDEX input.

  subroutine add_offp_index(this, offp_index)
    use integer_set_type_wavl
    class(index_map), intent(inout) :: this
    integer, intent(in) :: offp_index(:)
    type(integer_set) :: offp_set
    ASSERT(minval(offp_index) >= 1)
    ASSERT(maxval(offp_index) <= this%global_size)
    ASSERT(all((offp_index < this%first_gid) .or. (offp_index > this%last_gid)))
    call offp_set%add(offp_index) ! sort and remove duplicates
    call add_offp_index_set(this, offp_set)
  end subroutine

  !! Add off-process indices to an already initialized index map.
  !! This is an internal interface with trusted OFFP_SET input.

  subroutine add_offp_index_set(this, offp_set)

    use integer_set_type_wavl

    class(index_map), intent(inout) :: this
    type(integer_set), intent(inout) :: offp_set

    integer :: ierr, np, my_rank, rank, i, j, j1, n
    integer, allocatable :: last(:), offp_rank(:)
    integer, allocatable :: onp_count(:), offp_count(:)
    integer, allocatable :: onp_ranks(:), offp_ranks(:)

    !TODO? Allow extending an existing %offp_index
    INSIST(.not.allocated(this%offp_index))

    this%offp_index = offp_set
    call offp_set%clear
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    call MPI_Comm_size(this%comm, np, ierr)
    call MPI_Comm_rank(this%comm, my_rank, ierr)

    !! Determine which rank owns each off-process index (OFFP_RANK).
    !! OFFP_RANK will be ordered if OFFP_INDEX is ordered (we need this later).
    allocate(last(0:np-1), offp_rank(this%offp_size))
    call MPI_Allgather(this%last_gid, 1, MPI_INTEGER, last, 1, MPI_INTEGER, this%comm, ierr)
    rank = 0
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last(rank))
        rank = rank + 1
        ASSERT(rank < np)
      end do
      ASSERT(rank /= my_rank)
      offp_rank(j) = rank
    end do
    deallocate(last)

    !! Count the number of our off-process indices that are owned by each
    !! rank (OFFP_COUNT). This relies on the OFFP_RANK array being ordered.
    allocate(offp_count(0:np-1))
    j1 = 1
    do rank = 0, np-1
      do j = j1, size(offp_rank)
        if (offp_rank(j) > rank) exit
      end do
      offp_count(rank) = j - j1
      j1 = j
    end do
    deallocate(offp_rank)

    !! Distribute the number of our off-process indices that are owned by each
    !! rank to that rank. The result is the number of our on-process indices
    !! that are off-process on each rank (ONP_COUNT).
    allocate(onp_count(0:np-1))
    call MPI_Alltoall(offp_count, 1, MPI_INTEGER, onp_count, 1, MPI_INTEGER, this%comm, ierr)

    !! We now know which ranks need to communicate with which ranks in order
    !! to exchange data between off-process indices and their corresponding
    !! on-process indices. This all that is needed to define a virtual graph
    !! process topology.

    !! Compress the OFFP_COUNT array and generate the list of ranks that own
    !! our off-process indices: OFFP_RANKS is the list of incoming neighbors.
    n = count(offp_count > 0)
    allocate(offp_ranks(n), this%offp_counts(n))
    n = 0
    do rank = 0, np-1
      if (offp_count(rank) > 0) then
        n = n + 1
        offp_ranks(n) = rank
        this%offp_counts(n) = offp_count(rank)
      end if
    end do
    deallocate(offp_count)

    !! Compress the ONP_COUNT array and generate the list of ranks with off-
    !! process indices owned by us: ONP_RANKS is the list of outgoing neighbors.
    n = count(onp_count > 0)
    allocate(onp_ranks(n), this%onp_counts(n))
    n = 0
    do rank = 0, np-1
      if (onp_count(rank) > 0) then
        n = n + 1
        onp_ranks(n) = rank
        this%onp_counts(n) = onp_count(rank)
      end if
    end do
    deallocate(onp_count)

    !! Create the virtual topology. Here we use the OFFP_COUNTS and ONP_COUNTS
    !! as the edge weights (used if rank reordering is enabled). An alternative
    !! is to replace them with MPI_UNWEIGHTED.  Rank reordering is disabled.
    call MPI_Dist_graph_create_adjacent(this%comm, &
        size(offp_ranks), offp_ranks, this%offp_counts, &
        size(onp_ranks),  onp_ranks,  this%onp_counts, &
        MPI_INFO_NULL, .false., this%gather_comm, ierr)
    INSIST(ierr == MPI_SUCCESS)

    call MPI_Dist_graph_create_adjacent(this%comm, &
        size(onp_ranks),  onp_ranks,  this%onp_counts, &
        size(offp_ranks), offp_ranks, this%offp_counts, &
        MPI_INFO_NULL, .false., this%scatter_comm, ierr)
    INSIST(ierr == MPI_SUCCESS)

    deallocate(onp_ranks, offp_ranks)

    !! The components %OFFP_COUNTS, %OFFP_DISPLS, %ONP_COUNTS and %ONP_DIPSLS
    !! initialized here are meant for use with MPI_Neighbor_alltoallv. The
    !! component %ONP_INDEX will be used to fill the on-process buffer; the
    !! corresponding off-process buffer is off-process data array itself.

    !! Generate displacements into the on-process buffer for the start of the
    !! data for each neighbor rank.
    allocate(this%onp_displs, mold=this%onp_counts)
    if (size(this%onp_displs) > 0) then
      this%onp_displs(1) = 0
      do i = 2, size(this%onp_displs)
        this%onp_displs(i) = this%onp_displs(i-1) + this%onp_counts(i-1)
      end do
    end if

    !! Generate displacements into the off-process buffer for the start of
    !! the data for each neighbor rank.
    allocate(this%offp_displs, mold=this%offp_counts)
    if (size(this%offp_displs) > 0) then
      this%offp_displs(1) = 0
      do i = 2, size(this%offp_displs)
        this%offp_displs(i) = this%offp_displs(i-1) + this%offp_counts(i-1)
      end do
    end if

    !! Communicate the global off-process indices to their owning ranks.
    allocate(this%onp_index(sum(this%onp_counts)))
    call MPI_Neighbor_alltoallv(this%offp_index, this%offp_counts, this%offp_displs, MPI_INTEGER, &
        this%onp_index, this%onp_counts, this%onp_displs, MPI_INTEGER, this%scatter_comm, ierr)
    ASSERT(ierr == MPI_SUCCESS)
    this%onp_index = this%onp_index - this%first_gid + 1  ! map to local indices
    ASSERT(all(this%onp_index >= 1 .and. this%onp_index <= this%onp_size))

    ASSERT(gather_offp_verified(this))

  end subroutine add_offp_index_set

  !! Add %COUNTS and %DISPLS data used by gather/scatter methods.
  subroutine add_dist_coll_info(this)
    class(index_map), intent(inout) :: this
    integer :: n, ierr
    INSIST(.not.allocated(this%counts))
    n = merge(this%nproc, 0, this%is_root)
    allocate(this%counts(n), this%displs(n))
    call MPI_Gather(this%onp_size, 1, MPI_INTEGER, &
        this%counts, 1, MPI_INTEGER, this%root, this%comm, ierr)
    if (this%is_root) then
      this%displs(1) = 0
      do n = 2, this%nproc
        this%displs(n) = this%displs(n-1) + this%counts(n-1)
      end do
    end if
  end subroutine

  !! This function returns true if the gather_offp operation returns the
  !! expected result for integer data, and otherwise it returns false.
  !! Useful for testing with live index_map data.

  logical function gather_offp_verified(this) result(pass)
    type(index_map), intent(in) :: this
    integer :: j, ierr, onp_data(this%onp_size), offp_data(this%offp_size)
    do j = 1, this%onp_size
      onp_data(j) = global_index(this, j)
    end do
    call this%gather_offp(onp_data, offp_data)
    pass = all(offp_data == this%offp_index)
    call MPI_Allreduce(MPI_IN_PLACE, pass, 1, MPI_LOGICAL, MPI_LAND, this%comm, ierr)
  end function

  logical function defined(this)
    class(index_map), intent(in) :: this
    defined = .true.
  end function

end module index_map_type

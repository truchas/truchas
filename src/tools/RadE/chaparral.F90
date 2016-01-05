!!
!! CHAPARRAL
!!
!! Fortran 95 interface to the Chaparral vflib library (version 3.2)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chaparral

  implicit none
  public

  integer, parameter :: VF_MATRIX_EXCL_VIRT = 1
  integer, parameter :: VF_MATRIX_EXCL_DIAG = 2

  interface
    subroutine VF_Setup ()
    end subroutine

    subroutine VF_SetNumEnclosures (nenclosures)
      integer nenclosures
    end subroutine

    subroutine VF_SetMaxSurfaces (max_surfaces)
      integer max_surfaces
    end subroutine

    subroutine VF_CleanUp ()
    end subroutine

    integer function VF_DefineEnclosure (enclID, nonblocking, partial, asink, &
                                         npatches, global_ids, debug_level)
      integer, intent(in) :: nonblocking, partial, npatches, global_ids(*), debug_level
      character(len=*), intent(in) :: enclID
      double precision, intent(in) :: asink
    end function

    subroutine VF_RandomizeSurfacesOn ()
    end subroutine

    subroutine VF_RandomizeSurfacesOff ()
    end subroutine

    subroutine VF_DefineTopology (enclosure, geom_type, nfacets, &
         nnodes, x, y, z, c, f2p_map, nrotations, &
         x_mirror, y_mirror, z_mirror, bsp_depth, &
         bsp_length, spatial_tol, debug_level)
      integer, intent(in) :: enclosure, geom_type, nfacets, nnodes, c(*)
      integer, intent(in) :: f2p_map(*), nrotations, x_mirror, y_mirror, z_mirror
      integer, intent(in) :: bsp_depth, bsp_length, debug_level
      double precision, intent(in) :: spatial_tol, x(*), y(*), z(*)
    end subroutine

    subroutine VF_ResetTopology (enclosure)
      integer, intent(in) :: enclosure
    end subroutine

    subroutine VF_CalcHemicube (encl, hc_sub_divide, hc_resolution, hc_min_sep)
      integer, intent(in) :: encl, hc_sub_divide, hc_resolution
      double precision, intent(in) :: hc_min_sep
    end subroutine

    subroutine VF_CalcPairwise (encl, vis_nsamples, vis_sampling, &
                                mc_nsamples, mc_sampling, mc_tol1, mc_tol2)
      integer, intent(in) :: encl, vis_nsamples, vis_sampling, mc_nsamples, mc_sampling
      double precision, intent(in) :: mc_tol1, mc_tol2
    end subroutine

    subroutine VF_JitterOn ()
    end subroutine

    subroutine VF_JitterOff ()
    end subroutine

    subroutine VF_SmoothMatrix (encl, wt, tol, max_iter, symmetric, output)
      integer, intent(in) :: encl, max_iter, symmetric, output
      double precision, intent(in) :: wt, tol
    end subroutine

    subroutine VF_OutputMatrixSummaryBanner ()
    end subroutine

    subroutine VF_GetRowCounts (encl, mode, count)
      integer, intent(in)  :: encl, mode
      integer, intent(out) :: count(*)
    end subroutine

    subroutine VF_GetMatrix (encl, vf_cnt, vf_index, vf_data, vf_diag, vf_virt)
      integer, intent(in)  :: encl
      integer, intent(out) :: vf_cnt(*), vf_index(*)
      real,    intent(out) :: vf_data(*), vf_diag(*), vf_virt(*)
      optional vf_diag, vf_virt
    end subroutine
  end interface

end module chaparral

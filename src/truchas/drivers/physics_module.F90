!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!
!! NNC, 2 May 2012.  As part of removing old HT, the HEAT_CONDUCTION flag needed
!! to be moved somewhere.  Ultimately we want a "multi-process coordinator" that
!! is responsible for bundling together the various process kernels and managing
!! their interaction.  This flag belongs there; this is its _temporary_ home.
!!
!! Originally the HEAT_CONDUCTION flag was specific to the old solver, now it
!! refers to the new solver.
!!

module physics_module

  implicit none
  private

  logical, public, save :: heat_transport = .false.
  logical, public, save :: species_transport = .false.
  integer, public, save :: number_of_species = 0

end module physics_module

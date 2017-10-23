!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "f90_assert.fpp"

module face_output
   !----------------------------------------------------------------------------
   ! Purpose:
   !
   !   face output utilities
   !
   !   public interface:
   !
   !   call write_face_data
   !----------------------------------------------------------------------------
   use kinds, only: r8
   implicit none
   private
   
   integer, save, private :: unit
   ! public procedures
   integer, save, allocatable :: facelist(:)
   public :: write_face_data
   public :: open_face_file
   public :: close_face_file

CONTAINS

   !----------------------------------------------------------------------------
   subroutine write_face_data
      !-------------------------------------------------------------------------
      ! Purpose:
      !
      !    write face coordinates, temperatures, and fluxes to files, one per pe
      !-------------------------------------------------------------------------
   use diffusion_solver, only: ds_get_face_temp_copy,ds_get_face_flux_copy
   use unstr_mesh_type,  only: unstr_mesh 
   use mesh_manager,     only: unstr_mesh_ptr
   use time_step_module, only: t, dt, cycle_number 
     
   type(unstr_mesh),pointer :: mesh
   real(r8), allocatable::face_temp(:),face_flux(:) ! should these just be pointers?
   integer n,k,f
   real(r8) fx(3)
   character(1024) filename
   
   ! is this the best way to get the mesh pointer for this situation?
   mesh=>unstr_mesh_ptr('MAIN')

   
   
   ! get space for the copies of the face temperature and flux that are stored
   ! in the HT model
   allocate(face_temp(mesh%nface))
   allocate(face_flux(mesh%nface))

   ! get copies of the temperature and flux from the diffusion solver
   ! should these be replaced with functions that return pointers, rather than copies?
   call ds_get_face_temp_copy(face_temp) 
   call ds_get_face_flux_copy(face_flux)
   
      
   ! write the face data out. inlcude face number,cycle number, time, face coordinate,
   ! face temperature, face flux, and face normal (to interpret sign of flux)
   do n = 1,size(facelist)
    f = facelist(n)
    associate (fnode => mesh%fnode(mesh%xfnode(f):mesh%xfnode(f+1)-1))
      fx(1:3) = sum(mesh%x(:,fnode),dim=2) / size(fnode)
    end associate
    write(unit,'(2i10,100es20.10)') n,cycle_number,t,fx(:),face_temp(f),face_flux(f),mesh%normal(:,f)
   enddo 
   
   deallocate(face_temp)
   deallocate(face_flux)

   end subroutine write_face_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine open_face_file
   ! there are a couple goals in this routine (maybe it should be split into
   ! multiple routines, but oh well). First, it needs to make a list of
   ! all the faces that lie on the user specified bounding box.  Second, it
   ! needs to write out a file containing the zone temperatures for all zones
   ! inside the bounding box   
   use parallel_info_module, only : p_info
   use file_utility,     only: make_file_name
   use unstr_mesh_type,  only: unstr_mesh 
   use mesh_manager,     only: unstr_mesh_ptr
   use output_control, only: face_dump_bbox
   use time_step_module, only: t, dt, cycle_number 
   use zone_module, only: zone
    
   real(r8) bbox1(3),bbox2(3)
   type(unstr_mesh),pointer :: mesh
   character(1024) filename
   integer n,k,inside,numfaces,f,c
   real(r8) fx(3),cx(3),tol,mindist
   integer, allocatable::face_tag(:)
   integer count
   
   tol = 1e-6
    
   mesh=>unstr_mesh_ptr('MAIN')
    
   ! bbox1 is the x,y,z minimums of the bbox 
   bbox1(1) = face_dump_bbox(1)
   bbox1(2) = face_dump_bbox(2) 
   bbox1(3) = face_dump_bbox(3) 
   ! bbox2 is the x,y,z maximums of the bbox
   bbox2(1) = face_dump_bbox(4)
   bbox2(2) = face_dump_bbox(5) 
   bbox2(3) = face_dump_bbox(6) 
   
   ! max sure the user specified minimums are smaller than the maximums
   do k =1,3
    ASSERT(bbox2(k)>bbox1(k))
   enddo
   
   ! prepare to tag faces on the bbox
   allocate(face_tag(mesh%nface_onp)) 
   face_tag = 0
   numfaces = 0
   count = 0
   do n = 1,mesh%nface_onp
    associate (fnode => mesh%fnode(mesh%xfnode(n):mesh%xfnode(n+1)-1))
      fx(1:3) = sum(mesh%x(:,fnode),dim=2) / size(fnode)
    end associate
    ! if the face is in or on the bounding box, tag it
    ! assume it's inside
    inside = 1
    mindist = 1e50
    do k = 1,3
      mindist = min(abs(bbox1(k)-fx(k)),mindist)
      mindist = min(abs(bbox2(k)-fx(k)),mindist)
      ! if the face coordinate is outside the bbox, change its inside flag
      if (fx(k)<bbox1(k)-tol.or.fx(k)>bbox2(k)+tol) then
        inside = 0
      endif
    enddo 
    if (inside.eq.1) then
      count = count + 1
    endif
    ! if it's inside and on one of the boundaries, then it needs to be tagged
    if (inside.eq.1.and.mindist<tol) then
      face_tag(n) = 1
      numfaces = numfaces + 1
    endif
          
   enddo
   
   ! now actually create the list of faces that are on the bbox
   allocate(facelist(numfaces))
   numfaces = 0
   do n = 1,mesh%nface_onp
    if (face_tag(n).eq.1) then
      numfaces = numfaces+1
      facelist(numfaces) = n
    endif
   enddo
   deallocate(face_tag)


   ! write out the temperature data for all zones inside the bbox
   filename = make_file_name('tempdata',p_info%thisPE)

   open(newunit = unit, file = filename, status = 'replace')
   
   ! write out temperatures
   do c = 1,mesh%ncell_onp
      associate (cnode => mesh%cnode(mesh%xcnode(c):mesh%xcnode(c+1)-1))
        cx(1:3) = sum(mesh%x(:,cnode),dim=2) / size(cnode)
      end associate
      inside = 1
      do k =1,3
        if (cx(k)<bbox1(k).or.cx(k)>bbox2(k)) then
          inside = 0
        endif
      enddo
      if (inside.eq.1) write(unit,'(2i10,100es20.10)') c,cycle_number,t,cx(:),zone(c)%temp
   enddo 
   
   close(unit)



   ! open a file, one per pe, for the flux data to be written
   ! what would happen to this during a restart? Probably nothing good.
   
   filename = make_file_name('fluxdata',p_info%thisPE)

   open(newunit = unit, file = filename, status = 'replace')
   
   end subroutine open_face_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   subroutine close_face_file
    ! cleanup
    close(unit)
    deallocate(facelist)
   end subroutine close_face_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module face_output

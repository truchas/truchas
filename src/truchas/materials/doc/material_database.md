# The `material_database` Derived Type

The [`material_database_type`](../material_database_type.F90) module defines
the `material_database` derived type, which is a simple container for a
collection of `material` class objects. An instance is intended to serve
as an in-memory database of materials from which a client application can
select specific materials to include in a simulation. The
[`material_factory`](./material_factory.md) module provides a procedure for
populating a database object using a parameter list description of materials.

#### has_matl
```Fortran
type(material_database) :: db
db%has_matl(name)
```
Returns true if the database contains the material `name`; otherwise it
returns false.

#### matl_ref
```Fortran
type(material_database) :: db
db%matl_ref(name)
```
Returns a pointer to a `material` class object with the given `name`.
A null pointer is returned if the database does not contain the named
material.

#### add_matl
```Fortran
type(material_database) :: db
class(material), allocatable :: matl
call db%add_matl(matl)
```
Adds the material `matl` to the database. The allocation for `matl` is assumed
by the database object and `matl` is returned unallocated to the caller. It is
an error if the database already contains a material having the same name as
`matl`.

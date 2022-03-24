SPECIES_BC Namelist
====================

Overview
----------

The SPECIES_BC namelist is used to define boundary conditions for the species
diffusion model at external boundaries. Each instance of the namelist defines a
particular condition to impose on a species component over a subset of the
domain boundary. The boundary subset :math:`\Gamma` is specified using mesh face
sets. The namelist variable `face_set_ids`_ takes a list of face set IDs, and
the boundary condition is imposed on all faces belonging to those face sets.
Note that ExodusII mesh side sets are imported into Truchas as face sets with
the same IDs. The species component is specified using the `comp`_ namelist
variable.

:Required/Optional: Required
:Single/Multiple Instances: Multiple

The following types of boundary conditions can be defined. The outward unit
normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`.

- **Concentration**. A concentration Dirichlet condition for species component
  :math:`j`

  :math:`\phi_j = c` on :math:`\Gamma`

  is defined by setting type to ``'concentration'``. The boundary value
  :math:`c` is specified using either `conc`_ for a constant value, or
  `conc_func`_ for a function.

- **Flux**. A total concentration flux condition for species component :math:`j`

  :math:`-D_j\nabla\phi_j.\hat{n} = q` on :math:`\Gamma`,

  when the system does not include temperature as a dependent variable, or
  :math:`-D_j(\nabla\phi_j + S_j\nabla T).\hat{n} = q` on :math:`\Gamma` when it
  does, is defined by setting `type`_ to ``'flux'``. The concentration flux q is
  specified using either `flux`_ for a constant value, or `flux_func`_ for a
  function.

The specified species concentration boundary conditions are not allowed to
overlap, and they must completely cover the computational boundary.


Components
------------

.. contents::
   :local:


name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: A unique name used to identify a particular instance of this
              namelist.
:Type: string (31 characters max)
:Default: none


comp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The species component this boundary condition applies to.
:Type: integer
:Default: 1


face_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: A list of face set IDs that define the portion of the boundary
              where the boundary condition will be imposed.
:Type: integer list (32 max)
:Default: none


type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The type of boundary condition. The available options are:

              - ``'concentration'`` Concentration is prescribed on the boundary.
                Use `conc`_ or `conc_func`_ to specify
                its value.

              - ``'flux'`` Total outward concentration flux is prescribed on the
                boundary. Use `flux`_ or `flux_func`_
                to specify its value.

:Type: string
:Default: none


conc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The constant value of boundary concentration for a
              concentration-type boundary condition. To specify a function, use
              `conc_func`_ instead.
:Type: real
:Default: none


conc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The name of a :ref:`FUNCTION namelist
              <FUNCTION_Namelist/index:FUNCTION Namelist>` namelist defining a
              function that gives the boundary concentration for a
              concentration-type boundary condition. The function is expected to
              be a function of :math:`(t,x,y,z)`.
:Type: string
:Default: none


flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The constant value of the total outward boundary concentration
              flux for a flux-type boundary condition. To specify a function,
              use `flux_func`_ instead.
:Type: real
:Default: none


flux_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The name of a :ref:`FUNCTION namelist
              <FUNCTION_Namelist/index:FUNCTION Namelist>` namelist defining a
              function that gives the total outward boundary concentration flux
              for a flux-type boundary condition. The function is expected to be
              a function of :math:`(t,x,y,z)`.
:Type: string
:Default: none

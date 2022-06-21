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
normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`. The mass
flux is defined as :math:`\vec{f}_\mathrm{mass} \equiv -D_j(\nabla\phi_j +
S_j\nabla T)\cdot\hat{n}`.

External boundaries
^^^^^^^^^^^^^^^^^^^

- **Concentration**. A concentration Dirichlet condition for species component
  :math:`j`

  :math:`\phi_j = c` on :math:`\Gamma`

  is defined by setting type to ``'concentration'``. The boundary value
  :math:`c` is specified using either `conc`_ for a constant value, or
  `conc_func`_ for a function.

- **Flux**. A total concentration flux condition for species component :math:`j`

  :math:`\vec{f}_\mathrm{mass} = q` on :math:`\Gamma`,

  is defined by setting `type`_ to ``'flux'``. The concentration flux :math:`q`
  is specified using either `flux`_ for a constant value, or `flux_func`_ for a
  function.

- **Mass Transfer**: An external mass transfer flux condition on species
  component :math:`j`

  :math:`\vec{f}_\mathrm{mass} = k_j(\phi_j - \phi_{\infty,j})` on
  :math:`\Gamma`

  is defined by setting type to ``'mtc'``. The mass transfer coefficient
  :math:`k_j` for species :math:`j` is specified using either `mtc`_ or
  `mtc_func`_. The ambient concentration :math:`\phi_{\infty,j}` for species
  :math:`j` is specified using either `ambient_conc`_ or `ambient_conc_func`_.

The specified species concentration boundary conditions are not allowed to
overlap, and they must completely cover the computational boundary.

Internal interfaces
^^^^^^^^^^^^^^^^^^^
Internal interfaces are merely coincident pairs of conforming external mesh
boundaries. These are modifications to the mesh created by Truchas and are
defined using the :ref:`Interface_Side_Sets<M_ISS>` parameter from the
:ref:`MESH<MESH_Namelist>` namelist. Only the side set IDs referenced there
can be used in the definition of the following interface condition.

- **Interface Mass Transfer**: An interface mass transfer condition models
  mass transfer across an imperfect contact between two bodies or across a
  thin subscale material layer lying along an interface :math:`\Gamma`.
  It imposes continuity of the mass flux :math:`\vec{f}_\mathrm{mass}`
  across the interface :math:`\Gamma` and gives this flux as

      :math:`\vec{f}_\mathrm{mass} = k_j[\phi_j]` on :math:`\Gamma`
  
  where :math:`[\phi_j]` is the jump in :math:`\phi_j` across :math:`\Gamma`
  in the direction :math:`\hat{n}`. It is defined by setting `type`_ to
  ``'interface-mtc'``. The mass transfer coefficient :math:`k_j` for species
  :math:`j` is specified using either `mtc`_ for a constant value, or
  `mtc_func`_ for a function.

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

              - ``'concentration'``: Concentration is prescribed on the
                boundary. Use `conc`_ or `conc_func`_ to specify its value.

              - ``'flux'``: Total outward concentration flux is prescribed on
                the boundary. Use `flux`_ or `flux_func`_ to specify its value.

              - ``'mtc'``: External mass transfer condition. Use `mtc`_ or
                `mtc_func`_ to set the mass transfer coefficient, and
                `ambient_conc`_ or `ambient_conc_func`_ to set the ambinet
                concentration.

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


mtc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The constant value of the mass transfer coefficient for a mass
              transfer-type boundary condition. To specify a function, use
              `mtc_func`_ instead.
:Type: real
:Default: none


mtc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The name of a :ref:`FUNCTION namelist
              <FUNCTION_Namelist/index:FUNCTION Namelist>` namelist defining a
              function that gives the mass transfer coefficient for a mass
              transfer-type boundary condition. The function is expected to be a
              function of :math:`(t,x,y,z)` for an external condition, and a
              function of :math:`(\phi_j,t,x,y,z)` for an interface condition.
              In the latter case :math:`\phi_j` is taken to be the maximum of
              :math:`\phi_j` on either side of the interface.
:Type: string
:Default: none


ambient_conc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The constant value of the ambient concentration for species
              component :math:`j` for a mass transfer-type boundary condition.
              To specify a function, use `ambient_conc_func`_ instead.
:Type: real
:Default: none


ambient_conc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Description: The name of a :ref:`FUNCTION namelist
              <FUNCTION_Namelist/index:FUNCTION Namelist>` namelist defining a
              function that gives the ambient concentration for species
              component :math:`j` for a mass transfer-type boundary condition.
              The function is expected to be a function of
              :math:`(t,x,y,z)`.
:Type: string
:Default: none

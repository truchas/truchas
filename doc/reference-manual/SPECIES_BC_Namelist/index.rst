.. _SPECIES_BC_Namelist:

.. toctree::
   :maxdepth: 1

SPECIES_BC Namelist
====================

Overview
----------
The SPECIES_BC namelist is used to define boundary conditions for the species diffusion model at external boundaries. Each instance of the namelist defines a particular condition to impose on a species component over a subset of the domain boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets. The namelist variable :ref:`face_set_ids<S_BC_FSI>` takes a list of face set IDs, and the boundary condition is imposed on all faces belonging to those face sets. Note that ExodusII mesh side sets are imported into Truchas as face sets with the same IDs. The species component is specified using the :ref:`comp<S_BC_comp>` namelist variable.

The following types of boundary conditions can be defined. The outward unit normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`.

* **Concentration**. A concentration Dirichlet condition for species component :math:`j`

:math:`\phi_j = c` on :math:`\Gamma`

is defined by setting type to **"concentration"**. The boundary value :math:`c` is specified using either :ref:`conc<S_BC_C>` for a constant value, or :ref:`conc_func<S_BC_CF>` for a function.

* **Flux**. A total concentration flux condition for species component :math:`j`

:math:`-D_j\nabla\phi_j.\hat{n} = q` on :math:`\Gamma`,

when the system does not include temperature as a dependent variable, or 
:math:`-D_j(\nabla\phi_j + S_j\nabla T).\hat{n} = q` on :math:`\Gamma`
when it does, is defined by setting :ref:`type<S_BC_T>` to **"flux"**. The concentration flux q is specified using either :ref:`flux<S_BC_F>` for a constant value, or :ref:`flux_func<S_BC_FF>` for a function.

The specified species concentration boundary conditions are not allowed to overlap, and they must completely cover the computational boundary.

SPECIES_BC Namelist Features
----------------------------
| **Required/Optional        :** Required 
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`name<S_BC_name>`
* :ref:`face_set_ids<S_BC_FSI>`
* :ref:`comp<S_BC_comp>`
* :ref:`type<S_BC_T>`
* :ref:`conc<S_BC_C>`
* :ref:`conc_func<S_BC_CF>`
* :ref:`flux<S_BC_F>`
* :ref:`flux_func<S_BC_FF>`

.. _S_BC_name:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist.
| **Type**        : string (31 characters max)
| **Default**     : none

.. _S_BC_comp:

comp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The species component this boundary condition applies to.
| **Type**        : integer
| **Default**     : 1

.. _S_BC_FSI:

face_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face set IDs that define the portion of the boundary where the boundary condition will be imposed.
| **Type**        : integer list (32 max)
| **Default**     : none

.. _S_BC_T:

type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of boundary condition. The available options are:

* **concentration** Concentration is prescribed on the boundary. Use :ref:`conc<S_BC_C>` or :ref:`conc_func<S_BC_CF>` to specify its value.
* **flux** Total outward concentration flux is prescribed on the boundary. Use :ref:`flux<S_BC_F>` or :ref:`flux_func<S_BC_FF>` to specify its value.

| **Type**        : string
| **Default**     : none

.. _S_BC_C:

conc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of boundary concentration for a concentration-type boundary condition. To specify a function, use :ref:`conc_func<S_BC_CF>` instead. 
| **Type**        : real
| **Default**     : none

.. _S_BC_CF:

conc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the boundary concentration for a concentration-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`. 
| **Type**        : string
| **Default**     : none

.. _S_BC_F:

flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the total outward boundary concentration flux for a flux-type boundary condition. To specify a function, use :ref:`flux_func<S_BC_FF>` instead. 
| **Type**        : real
| **Default**     : none

.. _S_BC_FF:

flux_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the total outward boundary concentration flux for a flux-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`. 
| **Type**        : string
| **Default**     : none




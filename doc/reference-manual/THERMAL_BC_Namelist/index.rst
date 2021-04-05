.. _THERMAL_BC_Namelist:

.. toctree::
   :maxdepth: 1

THERMAL_BC Namelist
====================

Overview
------------

The THERMAL_BC namelist is used to define boundary conditions for the heat transfer model at external boundaries and internal interfaces. Each instance of the namelist defines a particular condition to impose over a subset of the domain boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets. The namelist variable :ref:`face_set_ids<TB_FSI>` takes a list of face set IDs, and the boundary condition is imposed on all faces belonging to those face sets. Note that ExodusII mesh side sets are imported into Truchas as face sets with the same IDs.

External boundaries
^^^^^^^^^^^^^^^^^^^^^^
The following types of external boundary conditions can be defined. The outward unit normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`.

.. _boundary_conditions_option:
.. csv-table:: Boundary Conditions 
   :header: "Option", "Equation", "Description"
   :class: tight-table
   :widths: 1 3 4

   "Temperature", ":math:`T = T_b` on :math:`\Gamma`", "A temperature Dirichlet condition is defined by setting :ref:`type<TB_Type>` to **temperature**. The boundary value :math:`T_b` is specified using either :ref:`temp<TB_Temp>` for a constant value, or :ref:`temp_func<TB_TempFunc>` for a function."
   "Total Flux", "| :math:`\kappa \Delta T.\hat{n}`
   | :math:`= q_b` on :math:`\Gamma`","A heat flux condition is defined by setting :ref:`type<TB_Type>` to **flux**. The heat flux :math:`{q_b}` is specified using either :ref:`flux<TB_Flux>` for a constant value, or :ref:`flux_func<TB_FluxFunc>` for a function."
   "Heat Transfer","| :math:`-\kappa\Delta T.\hat{n}`
   | :math:`= \alpha(T - T_{\infty})` on :math:`\Gamma`","An external heat transfer flux condition is defined by setting :ref:`type<TB_Type>` to **htc**. The heat transfer coefficient :math:`\alpha` is specified using either :ref:`htc<TB_htc>` for a constant value, or :ref:`htc_func<TB_htcFunc>` for a function, and the ambient temperature :math:`T_{\infty}` is specified using either :ref:`ambient_temp<TB_AT>` for a constant value, or :ref:`ambient_temp_func<TB_ATF>` for a function."
   "Ambient Radiation", "| :math:`-\kappa\Delta T.\hat{n}` :math:`= \epsilon\sigma((T-T_o)^4 - (T_{\infty} - T_o)^4)` on :math:`\Gamma`", "A simple ambient thermal radiation condition is defined by setting :ref:`type<TB_Type>` to **radiation**. The emissivityis specified using either :ref:`emissivity<TB_emis>` for a constant value or :ref:`emissivity_func<TB_emisFunc>` for a function, and the temperature of the ambient environment :math:`T_{\infty}` is specified using either :ref:`ambient_temp<TB_AT>` for a constant value, or :ref:`ambient_temp_func<TB_ATF>` for a function. Here :math:`\sigma` is the Stefan-Boltzmann constant and :math:`T_0` is the absolute-zero temperature, both of which can be redefined if the problem units differ from the default SI units using the :ref:`Stefan_Boltzmann<PhyCo_SB>` and :ref:`Absolute_Zero<PhyCo_AZ>` components of the :ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist."

The specified boundary conditions are not generally allowed to overlap. It is not permitted, for example, to imposed both a temperature and a flux condition on the same part of boundary. The one exception is that heat transfer and ambient radiation conditions can be superimposed; the net flux in this case will be the sum of the heat transfer and radiation fluxes.

It is also generally required that the specified boundary conditions completely cover the computational boundary. However, when enclosure radiation systems are present no boundary condition need be imposed on any part of the boundary that belongs to an enclosure. Either temperature or heat transfer conditions may still be imposed there, and in the latter case the net heat flux is the sum of the heat transfer and radiative (from enclosure radiation) fluxes.

Internal interfaces
^^^^^^^^^^^^^^^^^^^^^
Internal interfaces are merely coincident pairs of conforming external mesh boundaries. These are modifications to the mesh created by Truchas and are defined using the :ref:`Interface_Side_Sets<M_ISS>` parameter from the :ref:`MESH<MESH_Namelist>` namelist. Only the face set IDs referenced there can be used in the definition of the following interface conditions. The following types of internal interface conditions can be defined.

* **Interface Heat Transfer** An interface heat transfer condition models heat transfer across an imperfect contact between two bodies or across a thin subscale material layer lying along an interface :math:`\Gamma`. It imposes continuity of the heat flux :math:`−\kappa\Delta T·\hat{n}` across the interface :math:`\Gamma` and gives this flux as

:math:`-\kappa\Delta T.\hat{n} = -\alpha[T]` on :math:`\Gamma`,

where **[T]** is the jump in **T** across :math:`\Gamma` in the direction :math:`\hat{n}`. It is defined by setting :ref:`type<TB_Type>` to **interface-htc**. The heat transfer co-efficient :math:`\alpha` is specified using either :ref:`htc<TB_htc>` for a constant value, or :ref:`htc_func<TB_htcFunc>` for a function.

* **Gap Radiation**. A gap radiation condition models radiative heat transfer across a thin open gap lying along an interface :math:`\Gamma`. It imposes continuity of the heat flux :math:`-\kappa\Delta T.\hat{n}` across :math:`\Gamma` and gives the flux as

:math:`-\kappa\Delta T.\hat{n} = \epsilon_{\Gamma}\sigma((T_{-}-T_o)^4-(T_{+}-T_o)^4)` on :math:`\Gamma`,

where :math:`T_{-}` and :math:`T_{+}` denote the values of :math:`T` on the inside and outside gap surfaces with respect to the normal :math:`\hat{n}` to :math:`\Gamma`. It is defined by setting :ref:`type<TB_Type>` to **gap_radiation**. The gap emissivity :math:`\epsilon_{\Gamma}` is specified using either :ref:`emissivity<TB_emis>` for a constant value, or :ref:`emissivity_func<TB_emisFunc>` fir a function. The effective gap emissivity :math:`\epsilon_{\Gamma}` depends on the emissivities :math:`\epsilon_{-}` and :math:`\epsilon_{+}` of the surfaces on either side of the gap and is given by

:math:`\epsilon_{\Gamma} = \frac{\epsilon_{-}\epsilon_{+}}{\epsilon_{-} + \epsilon_{+} - \epsilon_{-}\epsilon_{+}}`

The value of the Stefan-Boltzmann constant :math:`\sigma` and the absolute-zero temperature :math:`T_0` can be redefined if the problem units differ from the default SI units using the :ref:`Stefan_Boltzmann<PhyCo_SB>` and :ref:`Absolute_Zero<PhyCo_AZ>` components of the :ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist.

THERMAL_BC Namelist Features
----------------------------
| **Required/Optional        :** Required
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`name<TB_N>`
* :ref:`face_set_IDs<TB_FSI>`
* :ref:`type<TB_type>`
* :ref:`temp<TB_temp>`
* :ref:`temp_func<TB_tempFunc>`
* :ref:`flux<TB_flux>`
* :ref:`flux_func<TB_fluxFunc>`
* :ref:`htc<TB_htc>`
* :ref:`htc_func<TB_htcFunc>`
* :ref:`ambient_temp<TB_AT>`
* :ref:`ambient_temp_func<TB_ATF>`
* :ref:`emissivity<TB_emis>`
* :ref:`emissivity_func<TB_emisFunc>`

.. _TB_N:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist.
| **Type**        : string (31 characters max)
| **Default**     : none

.. _TB_FSI:

face_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face set IDs that define the portion of the boundary where the boundary condition will be imposed.
| **Type**        : integer list (32 max)
| **Default**     : none

.. _TB_Type:

type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of boundary condition. The available options are:

.. csv-table:: Type of Boundary Conditions 
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 3 

   "**temperature**", "Temperature is prescribed on the boundary. Use :ref:`temp<TB_temp>` or :ref:`temp_func<TB_tempFunc>` to specify its value."
   "**flux**", "Outward heat flux is prescribed on the boundary. Use :ref:`flux<TB_flux>` or :ref:`flux_func<TB_fluxFunc>` to set its value."
   "**htc**","External heat transfer condition. Use :ref:`htc<TB_htc>` or :ref:`htc_func<TB_htcFunc>` to set the heat transfer coefficient, and :ref:`ambient_temp<TB_AT>` or :ref:`ambient_temp_func<TB_ATF>` to set the ambient temperature."
   "**radiation**","A simple ambient thermal radiation condition. Use :ref:`emissivity<TB_emis>` or :ref:`emissivity_func<TB_emisFunc>` to set the emissivity, and :ref:`ambient_temp<TB_AT>` or :ref:`ambient_temp_func<TB_ATF>` to set the temperature of the ambient environment."
   "**interface-htc**", "An internal interface heat transfer condition. Use :ref:`htc<TB_htc>` or :ref:`htc_func<TB_htcFunc>` to set the heat transfer coefficient."
   "**gap-radiation**","A gap thermal radiation condition. Use :ref:`emissivity<TB_emis>` or :ref:`emissivity_func<TB_emisFunc>` to set the emissivity."

| **Type**        : string
| **Default**     : none

.. _TB_temp:

temp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of boundary temperature for a temperature-type boundary condition. To specify a function, use :ref:`temp_func<TB_tempFunc>` instead.
| **Type**        : real
| **Default**     : none

.. _TB_tempFunc:

temp_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the boundary temperature for a temperature-type boundary condition. The function is expected to be a function of(t,x,y,z).
| **Type**        : string
| **Default**     : none

.. _TB_flux:

flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the outward boundary heat flux for a flux-type boundary condition. To specify a function, use :ref:`flux_func<TB_fluxFunc>` instead.
| **Type**        : real
| **Default**     : none

.. _TB_fluxFunc:

flux_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the outward boundary heatflux for a flux-type boundary condition. The function is expected to be a function of (t,x,y,z).
| **Type**        : string
| **Default**     : none

.. _TB_htc:

htc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the heat transfer coefficient for either an external or interface heat transfer-type boundary condition. To specify a function, use :ref:`htc_func<TB_htcFunc>` instead.
| **Type**        : real
| **Default**     : none

.. _TB_htcFunc:

htc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the heat transfer coefficient for either an external or interface heat transfer-type boundary condition. The function is expected to be a function of (t,x,y,z) for an external heat transfer-type boundary condition, and a function of (T,t,x,y,z) for an interface heat transfer-type boundary condition. In the latter case T is taken to be the maximum of the two temperatures on either side of the interface.
| **Type**        : string
| **Default**     : none

.. _TB_AT:

ambient_temp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the ambient temperature for external heat transfer or radiation-type boundary condition. To specify a function, use :ref:`ambient_temp_func<TB_ATF>` instead.
| **Type**        : real
| **Default**     : none

.. _TB_ATF:

ambient_temp_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the ambient temperature for external heat transfer or radiation-type boundary condition. The function is expected to be a function of (t,x,y,z).
| **Type**        : string
| **Default**     : none

.. _TB_emis:

emissivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of emissivity for a radiation-type boundary condition. To specify a function, use :ref:`ambient_temp_func<TB_ATF>` instead.
| **Type**        : real
| **Default**     : none

.. _TB_emisFunc:

emissivity_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the emissivity for a radiation-type boundary condition. The function is expected to be a function of (t,x,y,z).
| **Type**        : string
| **Default**     : none
.. _BODY_Namelist:

.. toctree::
   :maxdepth: 1

BODY Namelist
===============

Overview
------------
The BODY namelists define initial material distributions and conditions. The BODY namelists are processed in the order they appear, and identify the specified part of the computational domain not claimed by any preceding BODY namelist. Any “background” type BODY must be listed last. Each namelist is used to specify a geometry and initial state. The geometry is specified via the variables using an acceptable combination of :ref:`surface_name<B_SN>`, :ref:`axis<B_axis>`, :ref:`fill<B_fill>`, :ref:`height<B_height>`, :ref:`length<B_length>`, :ref:`mesh_material_number<B_MMN>`, :ref:`radius<B_radius>`, :ref:`rotation_angle<B_RA>`, :ref:`rotation_pt<B_RP>`, and :ref:`translation_pt<B_TP>`, hereafter referred to as geometry-type parameters. The initial state is specified using :ref:`material_number<B_MN>`, :ref:`velocity<B_velocity>`, :ref:`phi<B_phi>`, and :ref:`temperature<B_T>` or  :ref:`temperature_function<B_TF>`.

BODY Namelist Features
------------------------
| **Required/Optional        :** Required
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`axis<B_axis>`
* :ref:`fill<B_fill>`
* :ref:`height<B_height>`
* :ref:`length<B_length>`
* :ref:`material_name<B_MN>`
* :ref:`mesh_material_number<B_MMN>`
* :ref:`phi<B_phi>`
* :ref:`radius<B_radius>`
* :ref:`rotation_angle<B_RA>`
* :ref:`rotation_pt<B_RP>`
* :ref:`surface_name<B_SN>`
* :ref:`temperature<B_T>`
* :ref:`temperature_function<B_TF>`
* :ref:`translation_pt<B_TP>`
* :ref:`velocity<B_velocity>`


.. _B_axis:

axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The axis to be used for defining a cylinder or plane.
| **Type**        : string
| **Default**     : (none)
| **Valid Values**: `x`, `y`, `z` 

.. _B_fill:

fill
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The side of the surface to which material is to be inserted for this body.
| **Type**        : string
| **Default**     : `inside`
| **Valid Values**: `inside`, `outside`

.. _B_height:

height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Height of a cylinder body.
| **Physcial Dimension**: L
| **Type**        : real
| **Default**     : (none)
| **Valid Values**: (0.0, :math:`\infty`)

.. _B_length:

length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Length of each side of the box body, or the coefficients of an ellipse or ellipsoid body.
| **Physcial Dimension**: L
| **Type**        : real triplet
| **Default**     : (none)
| **Valid Values**: (0.0, :math:`\infty`)

.. _B_MN:

material_name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Name of the material, or material phase in the case of a multi-phase material, that occupies the volume of this body.
| **Type**        : string
| **Default**     : (none)


.. _B_MMN:

mesh_material_number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** :  List of material numbers (element block IDs) associated with the cells as defined in the mesh file. This parameter is only meaningful when surface_name = ’from mesh file’.
| **Type**        : integer list (16 max)
| **Default**     : (none)
| **Valid Values**: Existing material numbers in mesh file (if the mesh file is in Exodus/Genesis format, this is the mesh block number).

.. _B_phi:

phi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Initial value of the diffusion solver’s multi-component scalar field in the material body.
| **Physcial Dimension**: varies
| **Type**        : real vector of Num_Species (link broken, needs to be checked) values
| **Default**     : 0.0
| **Valid Values**: (-:math:`\infty`, :math:`\infty`)

.. _B_radius:

radius
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Radius of the geometric body (cylinder, sphere, ellipsoid).
| **Physcial Dimension**: L
| **Type**        : real 
| **Default**     : none
| **Valid Values**: (0.0, :math:`\infty`)

.. _B_RA:

rotation_angle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Angle (degrees) about the (x,y,z) axes this body is to be rotated. This variable is only supported for ’plane’ and ’cylinder’ body types.
| **Type**        : real triplet
| **Default**     : 0.0, 0.0, 0.0
| **Valid Values**: (-:math:`\infty`, :math:`\infty`)

.. _B_RP:

rotation_pt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Location of the point about which this body is to be rotated. This variable is only supported for ’plane’ and ’cylinder’ body types.
| **Physcial Dimension**: L
| **Type**        : real triplet
| **Default**     : 0.0, 0.0, 0.0
| **Valid Values**: (-:math:`\infty`, :math:`\infty`)

.. _B_SN:

surface_name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Type of surface characterizing the interface topology for this body. The available options are:

.. list-table::
   :widths: 10 50
   :header-rows: 1 
   :class: tight-table   
 
   * - Option
     - Description
   * - **background**
     -  background body will occupy all space which has not been claimed by previously listed BODY namelists. If provided, it must be the final BODY namelist provided. When specified, no other geometry-type parameters are relevant.
   * - **plane**
     - A plane is specified using :ref:`axis<B_axis>`, :ref:`rotation_angle<B_RA>`, :ref:`rotation_pt<B_RP>`, and :ref:`fill<B_fill>` to define the normal direction, and :ref:`translation_pt<B_TP>` to provide a point on the plane surface. The normal vector is an ’outward’ normal, such that the region defined is in the opposite direction of the normal vector unless :ref:`fill<B_fill>` = ’outside’. 
   * - **box**
     - A box is specified using :ref:`translation_pt<B_TP>` as the center, :ref:`length<B_length>` for the length of x, y, and z sides respectively, and :ref:`fill<B_fill>` to invert the shape. This shape does not support rotation. 
   * - **sphere**
     - A sphere is specified using :ref:`translation_pt<B_TP>` as the center, :ref:`radius<B_radius>`, and fill to invert the shape.  
   * - **ellipsoid** 
     - An ellipsoid of the form :math:`\frac{(x-x_o)^2}{{l_1}^2} + \frac{(y-y_o)^2}{{l_2}^2} + \frac{(z-z_o)^2}{{l_3}^2} <= 1` is specified using :ref:`translation_pt<B_TP>` as the center, :ref:`length <B_length>` for :math:`l_1, l_2, l_3`, and :ref:`fill<B_fill>` to invert the shape. This shape does not support rotation. 
   * - **ellipse** 
     - An infinitely long elliptic cylinder of the form :math:`\frac{(x-x_o)^2}{{l_1}^2} + \frac{(y-y_o)^2}{{l_2}^2} <= 1` is specified using :ref:`translation_pt<B_TP>` as the center, :ref:`length<B_length>` for :math:`l_1, l_2`, and :ref:`fill<B_fill>` to invert the shape. This shape does not support rotation, and will be aligned the z axis.
   * - **cylinder**
     - A cylinder is specified using :ref:`translation_pt<B_TP>` as the center of the base, :ref:`axis<B_axis>`, :ref:`rotation_angle<B_RA>`, and :ref:`rotation_pt<B_RP>` to define the orientation, :ref:`radius<B_radius>`, :ref:`height<B_height>`, and :ref:`fill<B_fill>` to invert the shape.


| **Type**        : string
| **Default**     : none


.. _B_T:

temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Initial constant temperature of the material body.
| **Physcial Dimension**: Θ
| **Type**        : real 
| **Default**     : none
| **Notes**: Either :ref:`temperature<B_T>` or :ref:`temperature_function<B_TF>` must specified, but not both.

.. _B_TF:

temperature_function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the initial temperature function for the material body. That function is expected to be a function of (x,y, z).
| **Type**        : string
| **Default**     : none
| **Notes**: Either :ref:`temperature<B_T>` or :ref:`temperature_function<B_TF>` must specified, but not both.

.. _B_TP:

translation_pt
^^^^^^^^^^^^^^^^^
| **Description** : Location to which each surface origin of this body is translated.
| **Physcial Dimension**: L
| **Type**        : real triplet
| **Default**     : 0.0, 0.0, 0.0
| **Valid Values**: (-:math:`\infty`, :math:`\infty`) 

.. _B_velocity:

velocity
^^^^^^^^^^^^^^^^^
| **Description** : Initial velocity of the material body.
| **Physcial Dimension**: L/T
| **Type**        : real triplet
| **Default**     : 0.0, 0.0, 0.0
| **Valid Values**: (-:math:`\infty`, :math:`\infty`) 


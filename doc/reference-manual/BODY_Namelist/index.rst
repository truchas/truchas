.. _BODY_Namelist:

.. toctree::
   :maxdepth: 1

BODY Namelist
===============
The BODY namelists define initial material distributions and initial values
of the field unknowns. They are processed in the order they appear in the
input file, and identify the specified part of the computational domain not
claimed by any preceding BODY namelist. Any `"background"`-type BODY must
appear last. Each namelist is used to specify a geometry and initial state.
The geometry is specified via the variables using an acceptable combination
of `surface_name`_, `axis`_, `fill`_, `height`_, `length`_,
`mesh_material_number`_, `radius`_, `rotation_angle`_, `rotation_pt`_, and
`translation_pt`_, hereafter referred to as geometry-type parameters. The
initial state is specified using `material_name`_, `velocity`_, `conc`_ or
`conc_func`_, and `temperature`_ or `temperature_function`_.


.. admonition:: Namelist Usage

   :Required/Optional: Required
   :Single/Multiple Instances: Multiple


Namelist Variables
------------------

.. contents::
   :local:


axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The axis to be used for defining a cylinder or plane.

:Type: string
:Default: none
:Valid Values: `"x"`, `"y"`, `"z"` 


fill
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The side of the surface to which material is to be inserted for this body.

:Type: string
:Default: `inside`
:Valid Values: `"inside"`, `"outside"`


height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Height of a cylinder body.

:Type: real
:Default: (none)
:Valid Values: :math:`{}>0`


length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Length of each side of the box body, or the coefficients of an ellipse or
ellipsoid body.

:Type: real triplet
:Default: (none)
:Valid Values: :math:`{}>0`


material_name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Name of the material, or material phase in the case of a multi-phase material,
that occupies the volume of this body.

:Type: string
:Default: none


mesh_material_number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
List of cell set IDs (Exodus element block IDs) that define the material body.
This variable is only meaningful when `surface_name`_ is `"from mesh file"`.

:Type: integer list
:Default: none


conc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The constant initial values of the multi-component scalar field in the material
body. To specify a function use `conc_func`_ instead. This is a vector-valued
variable with elements that correspond to the different components of the
scalar field. Each component must be assigned an initial value by either
this variable or `conc_func`_.

:Type: real
:Default: none


conc_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The names of :ref:`FUNCTION<FUNCTION_Namelist>` namelists defining functions
that give the initial values of the multi-component scalar field. The functions
are expected to be functions of :math:`(x,y,z)`. This is a vector-valued
variable with elements that correspond to the different components of the
scalar field. Each component must be assigned an initial value by either
this variable or `conc`_.

:Type: string
:Default: none


radius
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Radius of the geometric body (cylinder, sphere, ellipsoid).

:Type: real 
:Default: none
:Valid Values: :math:`{}>0`


rotation_angle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Angle (degrees) about the (x,y,z) axes this body is to be rotated. This
variable is only supported for `"plane:` and `"cylinder"` body types.

:Type: real triplet
:Default: :math:`0.0, 0.0, 0.0`


rotation_pt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Location of the point about which this body is to be rotated. This variable
is only supported for `"plane"` and `"cylinder"` body types.

:Type: real triplet
:Default: :math:`0.0, 0.0, 0.0`


surface_name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Type of surface characterizing the interface topology for this body.
The available options are:

.. list-table::
   :widths: 10 50
   :header-rows: 1 
   :class: tight-table   
 
   * - Value
     - Description
   * - **"background"**
     -  background body will occupy all space which has not been claimed by
        previously listed BODY namelists. If provided, it must be the final
        BODY namelist provided. When specified, no other geometry-type
        parameters are relevant.
   * - **"plane"**
     - A plane is specified using `axis`_, `rotation_angle`_, `rotation_pt`_,
       and `fill`_ to define the normal direction, and `translation_pt`_ to
       provide a point on the plane surface. The normal vector is an *outward*
       normal, such that the region defined is in the opposite direction of the
       normal vector unless `fill`_ is "outside". 
   * - **"box"**
     - A box is specified using `translation_pt`_ as the center, `length`_ for
       the length of x, y, and z sides respectively, and `fill`_ to invert the
       shape. This shape does not support rotation. 
   * - **"sphere"**
     - A sphere is specified using `translation_pt`_ as the center, `radius`_,
       and `fill`_ to invert the shape.  
   * - **"ellipsoid"** 
     - An ellipsoid of the form :math:`(x-x_0)^2/{l_1}^2 + (y-y_0)^2/{l_2}^2 +
       (z-z_0)^2/{l_3}^2 <= 1` is specified using `translation_pt`_ for the 
       center :math:`(x_0,y_0,z_0)`, `length`_ for :math:`l_1, l_2, l_3`, and
       `fill`_ to invert the shape. This shape does not support rotation. 
   * - **"ellipse"** 
     - An infinitely long elliptic cylinder of the form :math:`(x-x_0)^2/{l_1}^2
       + (y-y_0)^2/{l_2}^2 <= 1` is specified using `translation_pt` as the
       center, `length`_ for :math:`l_1, l_2`, and `fill`_ to invert the shape.
       This shape does not support rotation, and will be aligned the z axis.
   * - **"cylinder"**
     - A cylinder is specified using `translation_pt`_ for the center of the
       base, `axis`_, `rotation_angle`_, and `rotation_pt`_ to define the
       orientation, `radius`_, `height`_, and `fill`_ to invert the shape.

:Type: string
:Default: none


temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Initial constant temperature of the material body.

:Type: real 
:Default: none
:Notes: Either `temperature`_ or `temperature_function`_ must specified, but
        not both.


temperature_function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the
initial temperature function for the material body. That function is expected
to be a function of :math:`(x,y,z)`.

:Type: string
:Default: none
:Notes: Either `temperature`_ or `temperature_function`_ must specified, but
        not both.


translation_pt
^^^^^^^^^^^^^^^^^
Location to which each surface origin of this body is translated.

:Type: real triplet
:Default: :math:`0.0, 0.0, 0.0`


velocity
^^^^^^^^^^^^^^^^^
Initial velocity of the fluid material body.

:Type: real triplet
:Default: :math:`0.0, 0.0, 0.0`


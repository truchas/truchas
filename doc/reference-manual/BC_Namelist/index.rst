.. _BC_Namelist:

.. toctree::
   :maxdepth: 1

BC Namelist 
===========================

Overview
----------

The **BC** namelist is used to define boundary conditions for the solid mechanics model at external boundaries and internal material interfaces. The preferred method for specifying the mesh surface where a boundary condition applies is to reference a side set from the **ExodusII**-format mesh. Alternatively, the mesh surface can specified using a conic surface:

.. math::
 \color{blue} 0 = p(x,y,z) &= \color{blue}c_0 + c_x*{x} + c_y*y + c_z*z + c_{xx}*x^2 + c_{yy}*y^2 \\
 &\color{blue}+ c_{zz}*z^2 + c_{xy}*xy + c_{xz}*xz + c_{yz}*yz 
 :label: eq_3

A face belongs to the mesh surface whenever its centroid lies on this surface (see :ref:`Conic_Tolerance<BC_CT>`). The coefficients are specified using the **Conic_*** variables. Another method for solid mechanics is to specify nodes. The method is selected using :ref:`Surface_Name<BC_SN>`. The specified surface may also be restricted to lie within a bounding box.
 
BC Namelist Features
-----------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`BC_Name<BC_BCN>`
* :ref:`BC_Table<BC_BCT>` missing link in the latex manual
* :ref:`BC_Type<BC_BCType>`
* :ref:`BC_Value<BC_BCV>`
* :ref:`BC_Variable<BC_BCVar>`
* :ref:`Bounding_Box<BC_BB>`
* :ref:`Conic_Constant<BC_CC>`
* :ref:`Conic_Tolerance<BC_CT>`
* :ref:`Conic_X<BC_CX>`
* :ref:`Conic_XX<BC_CXX>`
* :ref:`Conic_XY<BC_CXY>`
* :ref:`Conic_XZ<BC_CXZ>`
* :ref:`Conic_Y<BC_CY>`
* :ref:`Conic_YY<BC_CYY>`
* :ref:`Conic_YZ<BC_CYZ>`
* :ref:`Conic_Z<BC_CZ>`
* :ref:`Conic_ZZ<BC_CZZ>`
* :ref:`Mesh_Surface<BC_MS>`
* :ref:`Node_Disp_Coords<BC_NDC>`
* :ref:`Surface_Name<BC_SN>`

.. _BC_BCN:

BC_Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A name used to identify a particular instance of this namelist.
| **Type**        : case-sensitive string
| **Default**     : none
| **Notes**       : This is optional and used for logging purposes only.

.. _BC_BCType:

BC_Type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of boundary condition.
| **Type**        : string
| **Default**     : Depends on :ref:`BC_Variable<BC_BCVar>`:
|   ’displacement’: ’x-traction’, ’y-traction’, ’z-traction’
| **Valid Values**: Depends on :ref:`BC_Variable<BC_BCVar>`:
|  ’displacement’:’x-traction’, ’y-traction’, ’z-traction’,’x-displacement’, ’y-displacement’, ’z-displacement’, ’normal-displacement’, ’normal-traction’, ’free-interface’, ’normal-constraint’, ’contact’
| **Notes**       : 
| * The solid mechanics displacement solution defaults to a traction-free surface with no displacement constraints (’x-traction’, ’y-traction’, and ’z-traction’ set to zero.)
| * The ’free-interface’, ’normal-constraint’ and ’contact’ types can only be specified for interfaces with gap elements.


.. _BC_BCV:

BC_Value
^^^^^^^^^^
| **Description** : Value(s) for the constant(s) used in this BC definition. See also :ref:`BC_Table<BC_BCT>`
| **Physcial Dimension** : Varies
| **Type**        : real (up to 24 values depending on :ref:`BC_Type<BC_BCType>`)
| **Default**     : 0.0 
| **Notes**       : The meaning of the items in the BC_Value list depends on the particular boundary condition.

.. csv-table:: 
   :header: "BC_Variable", "BC_Type", "Value Description", "Physcial Dimension", "Number of Values"
   :class: longtable
   :widths: auto

   "displacement", "x-displacement", "displacement", ":math:`L`", "1"
   "displacement", "y-displacement", "displacement", ":math:`L`", "1"
   "displacement", "z-displacement", "displacement", ":math:`L`", "1"
   "displacement", "x-traction", "traction (force/area)", ":math:`F/L^2`", "1"
   "displacement", "y-traction", "traction (force/area)", ":math:`F/L^2`", "1"
   "displacement", "z-traction", "traction (force/area)", ":math:`F/L^2`", "1"
   "displacement", "normal-displacement", "displacement", ":math:`L`", "1"
   "displacement", "free-interface", "not used", "-", "0"
   "displacement", "normal-constraint", "not used", "-", "0"
   "displacement", "contact", "not used", "-", "0"
   
.. _BC_BCVar:

BC_Variable
^^^^^^^^^^^^^
| **Description** : The name of the variable to which this boundary condition applies.
| **Type**        : case-insensitive string
| **Default**     : none
| **Valid Values**: "displacement"

.. _BC_BB:

Bounding_Box
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The extents in each dimension of a bounding box that restricts the extent of the mesh surface where the boundary condition is applied. This does not apply in the case :ref:`Surface_Name<BC_SN>` is **"node set"**.
| **Physical Dimension**: L
| **Type**        : A real array (:math:`x_{min}, x_{max}, y_{min}, y_{max}, z_{min}, z_{max}`)
| **Default**     : Unlimited in each dimension.

.. _BC_CC:

Conic_Constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_0` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0


.. _BC_CT:

Conic_Tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A mesh face is considered to lie on the conic surface when the absolute value of the conic polynomial (Equation :eq:`eq_3`) at the face centroid is less than the value of this parameter. Only relevant when using a conic polynomial to define the boundary condition surface.
| **Type**        : real 
| **Default**     : :math:`10^{-6}`
| **Valid Values**:  > 0
| **Notes**       : It is important to note that this is not a tolerance on the distance of a centroid from the conic surface, but merely a tolerance on the value of the conic polynomial. Its dimension depends on that of the coefficients in the polynomial. Most mesh generators place nodes on a bounding surface. For non-planar surfaces, this has the consequence that face centroids will not lie exactly on the surface, making the choice of this tolerance rather significant.

.. _BC_CX:

Conic_X
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_x` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CY:

Conic_Y
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_y` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CZ:

Conic_Z
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_z` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CXX:

Conic_XX
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{xx}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CYY:

Conic_YY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{yy}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CZZ:

Conic_ZZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{zz}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CXY:

Conic_XY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{xy}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CXZ:

Conic_XZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{xz}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_CYZ:

Conic_YZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the coefficient :math:`c_{yz}` in the conic polynomial (:eq:`eq_3`)
| **Type**        : real 
| **Default**     : 0.0

.. _BC_MS:

Mesh_Surface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Identifier of a side set defined in the **ExodusII**-format mesh. Only relevant when :ref:`Surface_Name<BC_SN>` is **"from mesh file"**.
| **Type**        : integer 
| **Default**     : none

.. _BC_NDC:

Node_Disp_Coords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : List of points that identify mesh nodes where a displacement boundary condition is applied. Up to 50 points can be specified as a list of (x,y,z) coordinates. Only relevant when :ref:`Surface_Name<BC_SN>` is **"node set"**.
| **Type**        : real
| **Physcial Dimension**     : :math:`L`

.. _BC_SN:

Surface_Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Selects the method of specifying the mesh surface where the boundary condition will be applied.
| **Type**        : case-insensitive string
| **Default**     : none
| **Valid Values**:

.. csv-table:: 
   :header: "Value", "Associated Variables"
   :class: tight-table
   :widths: auto

   "from mesh file", ":ref:`Mesh_Surface<BC_MS>`. Requires that the mesh is imported from an ExodusII-formatmesh file."
   "conic", ":ref:`Conic_Tolerance<BC_CT>`, :ref:`Conic_Constant<BC_CC>`, :ref:`Conic_X<BC_CX>`, :ref:`Conic_Y<BC_CY>`, :ref:`Conic_Z<BC_CZ>`, :ref:`Conic_XX<BC_CXX>`, :ref:`Conic_YY<BC_CYY>`, :ref:`Conic_ZZ<BC_CZZ>`, :ref:`Conic_XY<BC_CXY>`, :ref:`Conic_XZ<BC_CXZ>`, :ref:`Conic_YZ<BC_CYZ>`"
   "node set",":ref:`Node_Disp_Coords<BC_NDC>`"

.. _ENCLOSURE_RADIATION_Namelist:

.. toctree::
   :maxdepth: 1

ENCLOSURE_RADIATION Namelist 
==============================
  
Overview
------------

ENCLOSURE_RADIATION Namelist Features
-------------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple, one for each enclosure.

Components
------------
* :ref:`Name<ER_N>`
* :ref:`Enclosure_File<ER_EF>`
* :ref:`Coord_Scale_Factor<ER_CSF>`
* :ref:`Skip_Geometry_Check<ER_SGC>`
* :ref:`Ambient_Constant<ER_AC>`
* :ref:`Ambient_Function<ER_AF>`
* :ref:`Error_Tolerance<ER_ET>`
* :ref:`Precon_Method<ER_PM>`
* :ref:`Precon_Iter<ER_PI>`
* :ref:`Precon_Coupling_Method<ER_PCM>`
* :ref:`toolpath<ER_T>`

.. _ER_N:

Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : A unique name for this enclosure radiation system.
| **Type**        : string (31 characters max)
| **Default**     : none

.. _ER_EF:

Enclosure_File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The path to the enclosure file. This is interpreted relative to the Truchas input file unless this is an absolute path. If operating in moving enclosure mode (see :ref:`toolpath<ER_T>`) the file name is the base name for a collection of enclosure files.
| **Type**        : string (255 characters max)
| **Default**     : none
| **Valid Values**: (0, 0.1)
| **Notes**       : The genre program from the RadE tool suite can be used to generate this file.

.. _ER_CSF:

Coord_Scale_Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : An optional factor with which to scale the node coordinates of the enclosure surface.
| **Type**        : real
| **Default**     : 1.0
| **Valid Values**: :math:`\gt 0.0`
| **Notes**       : The faces of the enclosure surface must match faces from the Truchas mesh. If the coordinates of the mesh are being scaled, it is likely that the same scaling needs to be applied to the enclosure surface.

.. _ER_AC:

Ambient_Constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The constant value of the ambient environment temperature.
| **Physical dimension**: :math:`\Theta`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt` :ref:`Absolute_Zero<PhyCo_AZ>`
| **Notes**       : Either :ref:`Ambient_Constant<ER_AC>` or :ref:`Ambient_Function<ER_AF>` must be specified, but not both. Currently this is necessary even for full enclosures, although in that case the value will not be used and any value is acceptable.

.. _ER_AF:

Ambient_Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the ambient environment temperature function. That function is expected to be a function of :math:`t` alone.
| **Type**        : string
| **Default**     : none
| **Notes**       : Either :ref:`Ambient_Function<ER_AF>` or :ref:`Ambient_Constant<ER_AC>` must be specified, but not both. Currently this is necessary even for full enclosures, although in that case the value will not be used and any value is acceptable.

.. _ER_ET:

Error_Tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The error tolerance for the iterative solution of the linear radiosity system :math:`Aq=b`. Iteration stops when the approximate radiosity :math:`q` satisfies 

:math:`||b -Aq||_2 \lt \epsilon ||b||_2`

| **Type**        : real
| **Default**     : :math:`1.0e^{-3}`
| **Valid values**: :math:`\gt 0`
| **Notes**       : The Chebyshev iterative method is used when solving the radiosity system in isolation with given surface temperatures. However the usual case has the radiosity system as just one component of a larger nonlinear system that is solved by a Newton-like iteration, and this condition on the radiosity component is one necessary condition of the complete stopping criterion of the iteration.

.. _ER_T:

toolpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The name of a :ref:`TOOLPATH<TOOLPATH_Namelist>` namelist that defines a time-dependent motion that has been suitably partitioned. If this variable is specified it enables the moving enclosure mode of operation. This works in concert with the genre program.

.. _ER_PM:

Precon_Method (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Preconditioning method for the radiosity system.Use the default.
| **Type**        : string
| **Default**     : "jacobi"
| **Valid values**: "jacobi", "chebyshev"
| **Notes**       : The preconditioner for the fully-coupled heat transfer/enclosure radiation system NLK solver is built from smaller preconditioning pieces, one of which is a preconditioner for the radiosity system alone. The least costly and seemingly most effective is :math:`Jacobi`.

.. _ER_PI:

Precon_Iter (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The number of iterations of the :ref:`Precon_Method<ER_PM>` method to apply as the radiosity system preconditioner. Use the default.
| **Type**        : integer
| **Default**     : :math:`1`
| **Valid values**: :math:`\geq 1`

.. _ER_PCM:

Precon_Coupling_Method (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Method for coupling the radiosity and heat transfer system preconditioners. Use the default.
| **Type**        : string
| **Default**     : "backward GS"
| **Valid values**: "jacobi", "forward GS", "backward GS", "factorization"
| **Notes**       : There are several methods for combining the independent preconditionings of the radiosity system and heat transfer system to obtain a preconditioner for the fully-coupled system. If we view it as a block system with the the radiosity system coming first, the first three methods correspond to block Jacobi, forward block Gauss-Seidel, and backward Gauss-Seidel updates. The factorization method is an approximate Schur complement update, that looks like block forward Gauss-Seidel followed by the second half of block backward Gauss-Seidel.

.. _ER_SGC:

Skip_Geometry_Check (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Normally the geometry of the enclosure surface faces are compared with boundary faces of the heat conduction mesh to ensure they actually match. Setting this variable to false will disable this check, which may be necessary in some unusual use cases.
| **Type**        : logical
| **Default**     : false


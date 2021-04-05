.. _VFUNCTION_Namelist:

.. toctree::
   :maxdepth: 1

VFUNCTION Namelist
============================

Overview
------------
Similar to the :ref:`FUNCTION<FUNCTION_Namelist>` namelist, the :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist is used to define a vector-valued function that can be used in some situations where vector-valued data is needed, such as the specification of the boundary velocity in a flow boundary condition.

These are general vector-valued functions of one or more variables, :math:`y = (y_1,...,y_m) = f(x_1,...,x_n)`. The dimension of the value, the number of variables, and the unknowns they represent (i.e., time, position, temperature, etc.) all depend on the context in which the function is used, and that is described in the documentation of those namelists where these functions can be used. Currently only a single function type can be defined:

**Tabular Function.** This is a continuous, single-variable function :math:`y=f(s)` linearly interpolated from asequence of data points :math:`(s_i,y_i), i= 1,...,p`, with :math:`s_i \lt s_i+1`. The variable :math:`x_d` that is identified with :math:`s` is specified by :ref:`tabular_dim<VF_TDim>`.

VFUNCTION Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`Name<VF_Name>`
* :ref:`Type<VF_Type>`
* :ref:`Tabular_Data<VF_TDat>`
* :ref:`Tabular_Dim<VF_TDim>`

.. _VF_Name:

Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name by which this vector function can be referenced by other namelists.
| **Type**        : A case-sensitive string of up to 31 characters.
| **Default**     : None

.. _VF_Type:

Type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of function defined by the namelist.
| **Type**        : case-sensitive string 
| **Default**     : None
| **Valid values**: **tabular** for a tabular function

.. _VF_TDat:

Tabular_Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The table of values :math:`(s_i,y_i)` defining a tabular function :math:`y=f(s)`. Use :ref:`Tabular_Dim<VF_TDim>` to set the variable identified with :math:`s`.
| **Type**        : real array 
| **Default**     : None
| **Valid values**: This is a :math:`(m+ 1) × p` array that is most easily specified in the following manner

.. math::

   Tabular\_Data(:,1) =s_1, y_{11}, ..., y_{m1}
   
   Tabular\_Data(:,2) =s_2, y_{12}, ..., y_{m2}

      .

      .

      .

   Tabular\_Data(:,p) =s_p, y_{1p}, ..., y_{mp}

.. _VF_TDim:

Tabular_Dim
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The dimension in the :math:`m`-vector of independent variables that serves as the independent variable for the single-variable tabular function.
| **Type**        : integer
| **Default**     : 1


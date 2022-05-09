.. _VFUNCTION_Namelist:

VFUNCTION Namelist
============================
Similar to the :ref:`FUNCTION<FUNCTION_Namelist>` namelist, the VFUNCTION
namelist is used to define a vector-valued function that can be used in some
situations where vector-valued data is needed, such as the specification of
the boundary velocity in a flow boundary condition.

These are general vector-valued functions of one or more variables,
:math:`\mathbf{y} = (y_1,\ldots,y_m) = f(x_1,\ldots,x_n)`. The dimension of the value,
the number of variables, and the unknowns they represent (i.e., time,
position, temperature, etc.) all depend on the context in which the
function is used, and that is described in the documentation of those
namelists where these functions can be used.

The following function types are available.

* **Tabular Function.** This is a continuous, single-variable function
  :math:`y=f(s)` linearly interpolated from a sequence of data points
  :math:`(s_i,\mathbf{y}_i), i= 1,\ldots,p`, with :math:`s_i \lt s_{i+1}`.
  The variable :math:`x_d` that is identified with :math:`s` is specified
  by `tabular_dim`_.

* **Laser Irradiance Function.** This special-use type defines an alias to
  a laser irradiance function defined by a :ref:`TOOLHEAD<TOOLHEAD_Namelist>`
  namelist, which can be used by certain thermal boundary condition and
  source models.

.. note::

   :Required/Optional: Optional
   :Single/Multiple Instances: Multiple

Namelist Variables
--------------------------

.. contents::
   :local:


name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A unique name by which this vector function can be referenced by other
namelists.

:Type: case-sensitive string (31 characters max)
:Default: none


type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The type of function defined by the namelist.

:Type: case-sensitive string
:Default: none
:Valid values: "tabular", "toolhead-laser"


tabular_data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The table of p data points :math:`(s_i,\mathbf{y}_i)` defining a tabular function
:math:`\mathbf{y}=f(s)`. Use `tabular_dim`_ to set the independent variable
component identified with :math:`s`.

:Type: real array
:Default: none

.. note::
   This 1+m by p array is most easily specified point by point in this manner:

      | :fortran:`tabular_data(:,1) =` :math:`s_1, y_{11}, \ldots, y_{m1}`
      | :fortran:`tabular_data(:,2) =` :math:`s_2, y_{12}, \ldots, y_{m2}`
      | ...
      | :fortran:`tabular_data(:,p) =` :math:`s_p, y_{1p}, \ldots, y_{mp}`


tabular_dim
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The dimension in the :math:`n`-vector of independent variables that serves
as the independent variable for the single-variable tabular function.

:Type: integer
:Default: 1


toolhead
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The name of the :ref:`TOOLHEAD<TOOLHEAD_Namelist>` namelist that defines the
laser irradiance function that will be identified with this function. This is
a 3-vector valued function that depends on the variables :math:`(t, x, y, z)`.

:Type: case sensitive string
:Default: none

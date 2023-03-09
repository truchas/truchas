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

* **Shared Library Function** This is a function from a shared object library
  having a simple Fortran 77 or C compatible interface. Written in Fortran 77
  the function interface must look like

``subroutine myfun (v, p, r) bind(c)``

``double precision v(*), p(*), r(*)``

where **myfunc** can, of course, be any name. The equivalent C interface is

``void myfun (double v[], double p[], double r[]);``

The vector of variables :math:`v= (v_1,...,v_m)` is passed in the argument
:math:`v` and a vector of parameter values specified by `parameters`_ is
passed in the argument :math:`p`. The vector of variables
:math:`r=(r_1,...,r_d)` of dimension `dim`_ is returned. The path to the
library is given by `library_path`_ and the name of the function (myfun,
e.g.) is given by `library_symbol`_. Note that the bind(c) attribute on the
function declaration inhibits the Fortran compiler from mangling the function
name (by appending an underscore, for example) as it normally would.

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
:Valid values: "tabular", "toolhead-laser", "library"


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


library_path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The path to the shared object library that contains the function.

:Type: A string of up to 128 characters.
:Default: none


library_symbol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The symbol name of the function within the shared object file.

:Type: A string of up to 128 characters.
:Default: none
:Notes: Unless the Fortran function is declared with the BIND(C) attribute, which is the recommended practice, a Fortran compiler will almost always mangle the name of the function so that the symbol name is not quite the same as the name in the source code. Use the UNIX/Linux command-line utility **nm** to list the symbol names in the library file to determine the correct name to use here.

parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Optional parameter values to pass to the shared library function.

:Type: real vector of up to 16 values
:Default: None

dim
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Number of components of the return value of the shared library function.

:Type: integer
:Default: None

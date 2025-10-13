.. _COMPLEX_FUNCTION_Namelist:

COMPLEX_FUNCTION Namelist
=========================
Analogous to the :ref:`FUNCTION<FUNCTION_Namelist>` namelist, this namelist
is used to define a complex-valued function that can be used in some situations
where complex data is needed.

These are general complex-valued functions of one or more variables,
:math:`f:\mathbb{R}^m\to\mathbb{C}`. The number of variables, and the unknowns
they represent (e.g., time, position, etc.) depend on the context in which the
function is used, and that is described by the documentation of those namelists
where these functions can be used.

Currently the only type of function that is available is one from a shared
object library having a simple Fortran or C-compatible interface. Written
in Fortran the function interface must look like

.. code-block:: fortran

   subroutine myfun(v, p, r) bind(c)
   complex(kind(1.0d0)), intent(in)  :: v(*), p(*)
   complex(kind(1.0d0)), intent(out) :: r

where ``myfunc`` can, of course, be any name. The equivalent C interface is

.. code-block:: C

   #include <complex.h>
   void myfun(double v[], double p[], double complex *r);

The vector of variables :math:`v= (v_1,...,v_m)` is passed in the argument
``v`` and a vector of parameter values specified by `parameters`_ is passed
in the argument ``p``. The complex result is returned in the argument ``r``.
The path to the library is given by `library_path`_ and the name of the
function ("myfun", for example) is given by `library_symbol`_. Note that the
``bind(c)`` attribute on the function declaration inhibits the Fortran compiler
from mangling the function name (by appending an underscore, for example) as
it normally would.

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
:Valid values: "library"

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
:Notes: Unless the Fortran function is declared with the ``BIND(C)`` attribute,
   which is the recommended practice, a Fortran compiler will almost always
   mangle the name of the function so that the symbol name is not quite the
   same as the name in the source code. Use the UNIX/Linux command-line utility
   **nm** to list the symbol names in the library file to determine the correct
   name to use here.

parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Optional parameter values to pass to the shared library function.

:Type: real vector of up to 16 values
:Default: No parameters

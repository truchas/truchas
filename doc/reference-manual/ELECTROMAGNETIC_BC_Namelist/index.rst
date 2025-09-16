.. _ELECTROMAGNETIC_BC_Namelist:

.. toctree::
   :maxdepth: 2

ELECTROMAGNETIC_BC Namelist
===========================

The ELECTROMAGNETIC_BC namelist is used to define boundary conditions for the
electromagnetic (EM) models at external boundaries. Each instance of the
namelist defines a particular condition to impose over a subset of the domain
boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets.
The namelist variable `face_set_ids`_ takes a list of face set IDs, and the
boundary condition is applied at all faces belonging to those face sets. Note
that the mesh in view here is the one specified by the :ref:`EM_MESH
<EM_Mesh_Namelist>` namelist. Also note that ExodusII mesh side sets are
imported into Truchas as face sets with the same IDs.

.. admonishment:: Namelist Usage

   :Required/Optional: Required
   :Single/Multiple Instances: Multiple


Boundary Conditions
-------------------
The following types of boundary conditions can be defined. The outward unit
normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`.

Perfect electric conductor
^^^^^^^^^^^^^^^^^^^^^^^^^^
A perfect electric conductor condition, :math:`\hat{n}\times\vec{E}=0` on
:math:`\Gamma`, is defined by setting `type`_ to **"pec"**.

Perfect magnetic conductor
^^^^^^^^^^^^^^^^^^^^^^^^^^
A perfect magnetic conductor condition, :math:`\hat{n}\times\vec{H}=0` on
:math:`\Gamma`, is defined by setting `type`_ to **"pmc"**.

.. _em-bc-ih-field:

Magnetic induction source 
^^^^^^^^^^^^^^^^^^^^^^^^^
A tangential H-field condition :math:`\hat{n}\times\vec{H} =
\hat{n}\times\vec{H}_\text{ext}` on :math:`\Gamma`, where
:math:`\vec{H}_\text{ext}` is the external magnetic field defined by
the :ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist,
is defined by setting `type`_ to **"ih-hfield"**. This special purpose
boundary condition is used in induction heating simulations.

Robin
^^^^^^^^^^^^^^^^^^^^^^^^^
A general Robin condition

.. math::

   \hat{n}\times\nabla\times\vec{E} +
   \alpha\,\hat{n}\times\vec{E}\times\hat{n} =
   \hat{n}\times\vec{g}\times\hat{n}
   \text{ on $\Gamma$},
  
is defined by setting `type`_ to **"robin"**. The complex constant
:math:`\alpha` is specified by `alpha`_, and the complex function
:math:`\vec{g}(x)\in\mathbb{C}^3` is specified using either `g`_ for
a constant, or `g_func`_ for a function. This is only applicable to
the frequency domain EM model.

Impedance
^^^^^^^^^^^^^^^^^^^^^^^^^
This special form of a Robin condition approximates the effect of a
well-conducting bounding material, by choosing :math:`\alpha=
\sqrt{-i\mu_0\omega\sigma}` and :math:`\vec{g}=0`. The conductivity of the
bounding material :math:`\sigma` is specified by `sigma`_, and the angular
frequency :math:`\omega` is obtained from the configuration of the frequency
domain model. Set `type`_ to **"impedance"** to select this condition.


Waveguide port feed
^^^^^^^^^^^^^^^^^^^^^^^^^
This special form of a Robin condition models the entry and exit of the
fundamental TE\ :sub:`10` mode at the boundary of a truncated rectangular
waveguide. The geometry of the rectangular boundary is specified by variables
`center`_, `x_axis`_, `y_axis`_, `x_width`_, and `y_width`_, and the power of
the injected mode by `power`_. Set `type`_ to **"wg-port"** to select this
condition.

Namelist Variables
------------------ 
.. contents::
   :local:

name
^^^^^^^^^^^^^^^^^^
A unique name used to identify a particular instance of this namelist.

:Type: string
:Default: none

face_set_ids
^^^^^^^^^^^^^^^^^^
A list of face set IDs that define the portion of the boundary :math:`\Gamma`
where the boundary condition will be imposed.

:Type: integer list (32 max)
:Default: none

.. _em-bc-type:

type
^^^^^^^^^^^^^^^^^^
The type of boundary condition. Some types apply only to the frequency
domain (FD) model as noted. The available options are:

:Type: string
:Default: none
:Valid values:
.. list-table::
   :class: tight-table
   :header-rows: 1
   :widths: 1 5 5

   * - Value
     - Description
     - Associated variables
   * - **"pec"**
     - `Perfect electric conductor`_
     -
   * - **"pmc"**
     - `Perfect magnetic conductor`_
     -
   * - **"ih-hfield"**
     - `Magnetic induction source`_
     - :ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist
   * - **"robin"**
     - `Robin`_ (FD only)
     - `alpha`_, `g`_, `g_func`_
   * - **"impedance"**
     - `Impedance`_ (FD only)
     - `sigma`_
   * - **"wg-port"**
     - `Waveguide port feed`_ (FD only)
     - `center`_, `x_axis`_, `y_axis`_, `x_width`_, `y_width`_, `power`_,
       `phase`_


alpha
^^^^^^^^^^^^^^^^^^
The constant value of the :math:`\alpha` parameter in the Robin boundary
condition.

:Type: complex
:Default: none

g
^^^^^^^^^^^^^^^^^^
The constant value of the :math:`\vec{g}` parameter in the Robin boundary
condition. To specify a function use `g_func`_ instead.

:Type: complex 3-vector
:Default: none

g_func
^^^^^^^^^^^^^^^^^^
The name of a :ref:`COMPLEX_VFUNCTION<COMPLEX_VFUNCTION_Namelist>` namelist
that defines the :math:`\mathbb{C}^3`-valued function to use for the
:math:`\vec{g}`
parameter in the Robin boundary condition. The function is expected to be a
function of :math:`(x,y,z)`.

:Type: string
:Default: none

sigma
^^^^^^^^^^^^^^^^^^
The electric conductivity of the bounding material in the impedance boundary
condition.

:Type: real
:Default: none

.. _em-bc-power:

power
^^^^^^^^^^^^^^^^^^
The power of the input TE\ :sub:`10` mode in the waveguide port feed boundary
condition.

:Type: real
:Default: none
:Notes: No value needs to be specified for instances of the boundary condition
   that are managed by the :ref:`MICROWAVE_HEATING<MICROWAVE_HEATING_Namelist>`
   namelist, and if a value is specified it is ignored.

center
^^^^^^^^^^^^^^^^^^
The coordinate of the center of the boundary of the truncated rectangular
waveguide.

:Type: real 3-vector
:Default: none

x_axis
^^^^^^^^^^^^^^^^^^
The direction of the long edge of the rectangular boundary. Its length
is not significant.

:Type: real 3-vector
:Default: none

y_axis
^^^^^^^^^^^^^^^^^^
The direction of the short edge of the rectangular boundary. Its length
is not significant.

:Type: real 3-vector
:Default: none

x_width
^^^^^^^^^^^^^^^^^^
The length of the long edge of the rectangular boundary.

:Type: real
:Default: none

y_width
^^^^^^^^^^^^^^^^^^
The length of the short edge of the rectangular boundary.

:Type: real
:Default: none

phase
^^^^^^^^^^^^^^^^^^
The phase of the input TE\ :sub:`10` mode for the waveguide port feed boundary
condition. This is only truly significant when there are multiple waveguides
with differing phases.

:Type: real
:Default: 0.0


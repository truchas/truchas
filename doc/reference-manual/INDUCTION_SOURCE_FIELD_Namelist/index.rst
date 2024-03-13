.. _INDUCTION_SOURCE_FIELD_Namelist:

.. toctree::
   :maxdepth: 1

INDUCTION_SOURCE_FIELD Namelist 
===============================

The INDUCTION_SOURCE_FIELD namelist defines the external magnetic field that
serves as the applied source in the computation of Joule heat for induction
heating simulations. The field is a continuous wave (CW) function with
frequency :math:`f` of the form

.. math::
   \vec{H}_\text{ext}(\tau,\mathbf{x}) = \sin(2\pi f\tau)\,\vec{h}(\mathbf{x}),
   \quad\text{where}\quad
   \vec{h}(\mathbf{x}) = A\hat{e} + \sum_{i=1}^{m} \vec{h}_i(\mathbf{x})

is the superposition of a uniform field of strength :math:`A` oriented in
the unit direction :math:`\hat{e}`, and the magnetic fields :math:`\vec{h}_i`
produced by electric current in a set of induction coils. An
:ref:`INDUCTION_COIL <INDUCTION_COIL_Namelist>` namelist is used to specify
the parameters defining an induction coil.

Besides a steady CW source used for an entire induction heating simulation,
a CW source that changes abruptly at a sequence of heat transfer times
:math:`t_1,\ldots,t_n` may be defined by specifying addtional values for the
source parameters (:math:`f`, :math:`A`, and currents) to use at those times.

.. admonition:: Namelist Usage

   :Required/Optional: Required when :ref:`electromagnetics<PHYSICS_EM>` is enabled.
   :Single/Multiple Instances: Single

Namelist Variables
------------------

.. contents::
   :local:

orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The orientation :math:`\hat{e}` of the uniform field. For simplicity this
is currently restricted to one of the coordinates axes, :math:`\hat{x}`,
:math:`\hat{y}`, or :math:`\hat{z}`. This value also defines the orientation
of the induction coil axes.

:Type: string
:Default: `"z"`
:Valid Values: `"x"`, `"y"`, or `"z"`

.. _EM_ST:

times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional list of one or more heat transfer times :math:`t_1<\ldots<t_n`
when the CW source field changes.

:Type: real
:Default: none
:Valid Values: any strictly increasing sequence
:Notes: When this is specified, corresponding lists of values must also be
        specified for `frequency`_, `uniform_strength`_, and the
        :ref:`current<IC_Cur>` parameter in each :ref:`INDUCTION_COIL
        <INDUCTION_COIL_Namelist>` namelist. The number of values in those
        lists must be one more than the number of values specified here.

.. _EM_SF:

frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The frequency :math:`f` (cycles per unit time) of the sinusoidally-varying
source field. When `times`_ is specified, additional values to use at each
of those times must also be specified.

:Type: real
:Default: none
:Valid Values: any sequence of one or more postive values

.. _EM_US:

uniform_strength
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The strength :math:`A` of the uniform source field. The field is oriented
in the direction specified by `orientation`_. When `times`_ is specified,
additional values to use at each of those times must also be specified.

:Type: real
:Default: 0
:Valid Values: any sequence of one or more values
:Notes: A negative value essentially shifts the phase by 180 degrees.

        For reference, the magnitude of the magnetic field within an
        infinitely-long, finely-wound coil with linear current density
        :math:`I` is simply :math:`I`. The field magnitude at the center
        of a circular current loop of radius :math:`r` with current
        :math:`I` is :math:`I/{2r}`. In both cases the field is directed
        along the axis of the coil/loop.

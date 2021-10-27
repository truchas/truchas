.. _TOOLHEAD_Namelist:

TOOLHEAD Namelist
=================
The TOOLHEAD namelist is used to model an external energy source that is
attached to a moving toolhead like those employed in additive manufacturing
processes. Two choices of laser model are currently implemented. The motion
of the toolhead relative to the computational domain is specified using a
:ref:`TOOLPATH<TOOLPATH_Namelist>` namelist.

.. note::

   :Required/Optional: Optional
   :Single/Multiple Instances: Multiple

Namelist Variables
--------------------------

.. contents::
   :local:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A unique name used to identify a particular instance of this namelist.
Clients will reference the toolhead using this name.

:Type: case sensitive string (31 characters max)
:Default: none


toolpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The name of a :ref:`TOOLPATH<TOOLPATH_Namelist>` namelist that defines the
motion of the toolhead relative to the computational domain. The toolpath
position defines the reference point for the energy source associated with
the toolhead. The toolpath also specifies when the source is on and off via
the setting of flag 0: on (flag set), and off (flag clear).

:Type: case sensitive string (31 characters max)
:Default: none


laser_direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The direction of laser energy flow. The magnitude of this 3-vector is not
significant.

:Type: real 3-vector
:Valid values: any non-zero vector
:Default: none


laser_power
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The total power of the laser.

:Type: real
:Units: **E**/**T**
:Valid values: :math:`\ge0`
:Default: none


laser_time_const
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The exponential-decay time constant for the laser power employed when the
laser is turned on/off. Rather than the laser power changing discontinuously
at such a time :math:`t_0`, it fades in continuously to full power or out
to zero power with an exponential decay :math:`\exp((t_0-t)/\tau)` where
:math:`\tau` is the time constant of the decay. To have the power change
discontinuously, specify 0.

:Type: real
:Units: **T**
:Valid values: :math:`\ge0`
:Default: none

.. note::
   When the laser changes state between on/off this creates an abrupt change
   in forcing which results in a short-time transient which must be resolved
   by the time stepping procedure. Because of this, Truchas automatically
   reduces the step size by a significant factor whenever this occurs.
   In spite of this it can still pose significant numerical difficulties for
   time stepping. Judicious use of this exponential fade feature to soften
   the change will ease the time stepping without noticeably changing the
   long-time solution.

.. seealso::
   The :ref:`SIMULATION_CONTROL<SIMULATION_CONTROL_Namelist>` namelist
   provides a more brute force method of reducing the time step further,
   which *may* be helpful in some difficult cases.
   

laser_absorp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The laser absorption coefficient.  This is the fraction of the incident
power that is absorbed as heat.

:Type: real
:Valid values: :math:`[0.0,1.0]`
:Default: none


laser_type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the type of laser model used. See the following section for a
description of the available models and their associated namelist variables.

:Type: string
:Valid values: "gaussian" and "gaussian beam"
:Default: none

Laser Models
------------

Gaussian
^^^^^^^^
In the simple Gaussian model, the laser energy propagates along rays
parallel to the beam axis in the direction specified by `laser_direction`_.
The axis passes through the `toolpath`_ position, and the radiant flux
:math:`E_e` on any plane orthogonal to the axis is a Gaussian in the radius
:math:`r` from the axis:

.. math::

   E_e(r) \propto \frac{1}{2\pi\sigma^2}\exp\left(\frac{-r^2}{2\sigma^2}\right).

The model parameter :math:`\sigma` is specified using this namelist variable:

laser_sigma
+++++++++++
The value of :math:`\sigma`, which characterizes the size of the beam.
The relationship between :math:`\sigma` and the full-width-at-half-maximum
(FWHM) value of the beam is given by :math:`\text{FWHM} = 2\sqrt{2\ln2}\;\sigma
\approx 2.35\,\sigma`.

:Type: real
:Units: **L**
:Valid values: :math:`>0`
:Default: none

Gaussian Beam
^^^^^^^^^^^^^
The Gaussian beam model models a focused Gaussian laser beam whose energy
propagates along converging/diverging rays that are approximately parallel
to the beam axis with direction specified by `laser_direction`_. The focal
point of the beam lies on the axis and is located at the `toolpath`_ position.
The radiant flux :math:`E_e` on planes orthogonal to the axis is a Gaussian
in the radius :math:`r` from the axis, with a profile (amplitude and width)
that varies with distance :math:`z` of the plane from the focal point:

.. math::

   E_e(r,z) \propto \frac{2}{\pi w^2}\exp\left(\frac{-2r^2}{w^2}\right), \quad
   w(z) = w_0 \sqrt{1 + \left(\frac{z\lambda M^2}{\pi w_0^2}\right)^2}.

The model parameters :math:`w_0`, :math:`\lambda`, and :math:`M^2` are
specified using these namelist variables:

laser_waist_radius
++++++++++++++++++
The value of :math:`w_0`, which characterizes the size of the beam
at its focal point or waist. The relationship between :math:`w_0` and the
full-width-at-half-maximum (FWHM) value of the beam at its waist is given
by :math:`\text{FWHM} = \sqrt{2\ln2}\;w_0 \approx 1.18\,w_0`.

:Type: real
:Units: **L**
:Valid values: :math:`>0`
:Default: none

laser_wavelength
++++++++++++++++
The wavelength :math:`\lambda` of the laser radiation.

:Type: real
:Units: **L**
:Valid values: :math:`>0`
:Default: none

laser_beam_quality_factor
+++++++++++++++++++++++++
The unitless :dfn:`beam quality factor` :math:`M^2` which characterizes the
deviation of a laser beam from an ideal Gaussian beam. This is a widely-used
measure in the laser industry.

:Type: real
:Valid values: :math:`\ge1`
:Default: none

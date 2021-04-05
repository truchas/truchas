.. _VISCOPLASTIC_MODEL_Namelist:

.. toctree::
   :maxdepth: 1

VISCOPLASTIC_MODEL Namelist
============================

Overview
----------
The VISCOPLASTIC_MODEL namelist defines the viscoplastic model to be used for a particular solid material phase in material stress-strain calculations. Two viscoplastic models are available: a mechanical threshold stress (MTS) model and a power law model. When no model is given for a solid material phase, it is modeled as a purely elastic material. The models specify a relation for the effective plastic strain rate ̇as a function of temperature :math:`T` and von Mises stress :math:`\sigma`. For more details on the models see the Truchas Physics and Algorithms manual. Briefly:

**Power Law Model**. In the simple power law model, the strain rate relation is

.. math::
   :label: vpm_eq1
   
   \dot{\epsilon} = A*exp(-Q/RT)*\sigma^n

where A, n, Q, and R are parameters given by this namelist.

**MTS Model**. The MTS model uses the strain rate relation

.. math::
   :label: vpm_eq

   \dot{\epsilon} = \dot{\epsilon_{0i}} exp[-\frac{\mu b^3 g_{0i}}{kT}(1-(\frac{\mu_0}{\mu\sigma_i}(\sigma - \sigma_a))^{p_i})^{q_i}],  \mu = \mu_0 - \frac{D}{exp(T_0/T) - 1}

where :math:`\dot{\epsilon_{0i}}`, :math:`\dot{g_{0i}}`, b, k, D, :math:`\mu_0`, :math:`T_0`, :math:`\sigma_i`, :math:`\sigma_a`, :math:`p_i`, and :math:`q_i` are parameters given by this namelist. When :math:`\sigma \lt \sigma_a` we instead use :math:`\dot{\epsilon} = K\sigma^5`, where K is chosen to give continuity with the previous relation at :math:`\sigma = \sigma_a`. And when :math:`\sigma - \sigma_a \gt \mu\sigma_i/\mu_0` we take :math:`\dot{\epsilon} = \dot{\epsilon_{0i}}`.

VISCOPLASTIC_MODEL Namelist Features
----------------------------------------
| **Required/Optional        :** Optional; only relevant when :ref:`Solid_Mechanics<SOLID_MECHANICS_Namelist>` is true.
| **Single/Multiple Instances:** Multiple; at most one per solid material phase.

Components
------------
* :ref:`Phase <VPM_Phase>`
* :ref:`Model <VPM_Model>`
* :ref:`MTS_b <VPM_MTSb>`
* :ref:`MTS_d <VPM_MTSd>`
* :ref:`MTS_edot_0i <VPM_MTSe>`
* :ref:`MTS_g_0i <VPM_MTSg>`
* :ref:`MTS_k <VPM_MTSk>`
* :ref:`MTS_mu_0 <VPM_MTSm>`
* :ref:`MTS_p_i <VPM_MTSp>`
* :ref:`MTS_q_i <VPM_MTSq>`
* :ref:`MTS_sig_a <VPM_MTSsiga>`
* :ref:`MTS_sig_i <VPM_MTSsigi>`
* :ref:`MTS_temp_0 <VPM_MTSt>`
* :ref:`Pwr_Law_A <VPM_PLA>`
* :ref:`Pwr_Law_N <VPM_PLN>`
* :ref:`Pwr_Law_Q <VPM_PLQ>`
* :ref:`Pwr_Law_R <VPM_PLR>`

.. _VPM_Phase:

Phase
^^^^^^^^^^^^^^^^^^^

| **Description** : The name of the material :ref:`PHASE<MATERIAL_and_PHASE_Namelists>` to which this viscoplastic model applies.
| **Type**        : case-sensitive string
| **Default**     : none

.. _VPM_Model:

Model
^^^^^^^^^^^^^^^^^^^

| **Description** : The type of viscoplastic strain rate model.
| **Type**        : case-insensitive string
| **Default**     : none
| **Valid Values**: "MTS", "power law", "elastic"
| **Notes**       : The effect of the "elastic" option is equivalent to not specifying a viscoplastic model at all; it is provided as a convenience.

.. _VPM_MTSb:

MTS_b
^^^^^^^^^^^^^^^^^^^

| **Description** : Burgers vector length :math:`b` in :eq:`vpm_eq`
| **Physical Dimension**: :math:`L`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSd:

MTS_d
^^^^^^^^^^^^^^^^^^^

| **Description** : Constant :math:`D` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none

.. _VPM_MTSe:

MTS_edot_0i
^^^^^^^^^^^^^^^^^^^

| **Description** : Reference strain rate :math:`\dot{\epsilon_{0i}}` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`T^{-1}`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSg:

MTS_g_0i
^^^^^^^^^^^^^^^^^^^

| **Description** : Material constant :math:`g_{0i}` used in :eq:`vpm_eq`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSk:

MTS_k
^^^^^^^^^^^^^^^^^^^

| **Description** : Boltzmanns constant :math:`k` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`E/\Theta`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`
| **Note**        : Temperature should be expressed in Kelvin, or other temperature scale where 0 corresponds to absolute zero. If SI units are being used, :math:`k` should be :math:`1.38×10^{−23}`. Use a value appropriate to the units used in :eq:`vpm_eq`

.. _VPM_MTSm:

MTS_mu_0
^^^^^^^^^^^^^^^^^^^

| **Description** : Reference value :math:`\mu_0` for the temperature dependent shear modulus used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSp:

MTS_p_i
^^^^^^^^^^^^^^^^^^^

| **Description** : Exponent term :math:`p_i` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSq:

MTS_q_i
^^^^^^^^^^^^^^^^^^^

| **Description** : Exponent term :math:`q_i` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSsiga:

MTS_sig_a
^^^^^^^^^^^^^^^^^^^

| **Description** : The athermal stress term :math:`\sigma_a` used in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0`

.. _VPM_MTSsigi:

MTS_sig_i
^^^^^^^^^^^^^^^^^^^

| **Description** : A stress term :math:`\sigma_i` related to obstacles to dislocation motion in :eq:`vpm_eq`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_MTSt:

MTS_temp_0
^^^^^^^^^^^^^^^^^^^

| **Description** : Constant :math:`T_0` used in the temperature dependent shear modulus in :eq:`vpm_eq`
| **Physical Dimension**: :math:`\Theta`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_PLA:

Pwr_Law_A
^^^^^^^^^^^^^^^^^^^

| **Description** : Constant term :math:`A` in :eq:`vpm_eq1`
| **Physical Dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0`

.. _VPM_PLN:

Pwr_Law_N
^^^^^^^^^^^^^^^^^^^

| **Description** : Stress exponent term :math:`n` in :eq:`vpm_eq1`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _VPM_PLQ:

Pwr_Law_Q
^^^^^^^^^^^^^^^^^^^

| **Description** : Activation energy :math:`Q` in :eq:`vpm_eq1`
| **Physical Dimension**: :math:`E/mol`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0`

.. _VPM_PLR:

Pwr_Law_R
^^^^^^^^^^^^^^^^^^^

| **Description** : Gas constant :math:`R` in :eq:`vpm_eq1`
| **Physical Dimension**: :math:`E/(\Theta mol)`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`
| **Note**        : Temperature should be expressed in Kelvin, or other temperature scale where :math:`0` corresponds to absolute zero. Use the value for :math:`R` appropriate to the units used in :eq:`vpm_eq1`





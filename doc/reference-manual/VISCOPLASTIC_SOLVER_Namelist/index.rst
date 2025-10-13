VISCOPLASTIC_SOLVER Namelist
============================

The ``VISCOPLASTIC_SOLVER`` namelist sets parameters that are specific to the
viscoplasticity solver. This namelist is read whenever the :ref:`PHYSICS
Namelist<PHYSICS_Namelist/index:PHYSICS Namelist>` option
:ref:`solid_mechanics<PHYSICS_Namelist/index:solid_mechanics>` is enabled and a
:ref:`VISCOPLASTIC_MODEL
Namelist<VISCOPLASTIC_MODEL_Namelist/index:VISCOPLASTIC_MODEL Namelist>` exists.

:Required/Optional: Required when viscoplasticity is enabled (i.e., when a
   :ref:`VISCOPLASTIC_MODEL<VISCOPLASTIC_MODEL_Namelist/index:VISCOPLASTIC_MODEL Namelist>`
   namelist exists).
:Single/Multiple Instances: Single

.. contents:: Components
   :local:


strain_limit
^^^^^^^^^^^^^^^^^^^^^^^^^

This parameter controls the use of the ODE integrator in the plastic strain
calculation. When the plastic strain at an integration point is below this
limit, a single Heun step is taken. When the plastic strain is at or above this
limit, the solver requested by :ref:`solver
<VISCOPLASTIC_SOLVER_Namelist/index:solver (expert)>` is invoked with
an initial step size such that the predicted plastic strain delta does not
exceed the limit.

:Type: real
:Default: 1e-10
:Valid Values: :math:`\geq 0`

.. tip::

   This should be set to the minimum significant value of the plastic strain
   increment for a time step. If convergence seems poor when a viscoplastic
   material model is used, it may help to reduce this value.


abs_plastic_strain_tol
^^^^^^^^^^^^^^^^^^^^^^

The tolerance :math:`\epsilon` for the absolute error component of the plastic
strain error norm used by the nonlinear solver. If :math:`\delta u` is a plastic
strain field increment with reference plastic strain field :math:`u`,then this
error norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The relative error tolerance :math:`\eta` is given by `rel_plastic_strain_tol`_.

:Physical Dimension: :math:`\Theta`
:Type: real
:Default: 1e-12
:Valid Values: :math:`\geq 0`

.. note::
   The error norm is dimensionless and normalized.

.. note::
   For :math:`u_j` sufficiently small the norm approximates an absolute norm
   with tolerance :math:`\epsilon`, and for :math:`u_j` sufficiently large the
   norm approximates a relative norm with tolerance :math:`\eta`. If
   :math:`\epsilon = 0` then the norm is a pure relative norm and the
   plastic strain must be bounded away from 0.


rel_plastic_strain_tol
^^^^^^^^^^^^^^^^^^^^^^

The tolerance :math:`\eta` for the relative error component of the plastic
strain error norm used by the nonlinear solver. If :math:`\delta u` is a plastic
strain field increment with reference displacement field :math:`u`, then this
error norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The absolute error tolerance :math:`\epsilon` is given by `abs_plastic_strain_tol`_.

:Physcial Dimension: dimensionless
:Type: real
:Default: 1e-3
:Valid Values: (0, 1)

.. note::
   See the notes for `abs_plastic_strain_tol`_.


maximum_iterations
^^^^^^^^^^^^^^^^^^

Maximum allowed number of iterations of the nonlinear solver.

:Type: integer
:Default: 10
:Valid Values: :math:`[0,\infty)`


nlk_max_vectors
^^^^^^^^^^^^^^^

For the NLK method, the maximum number of acceleration vectors to be used.

:Type: integer
:Default: 3
:Valid Values: :math:`[0,\infty)`

rate_limit
^^^^^^^^^^^

This parameter controls the use of the ODE integrator in the plastic strain
calculation. When the relative rate of change of the plastic strain at an
integration point is below this limit, a single Heun step is taken. Otherwise,
(and if the `strain_limit`_ is met), the solver requested by
:ref:`solver<VISCOPLASTIC_SOLVER_Namelist/index:solver (expert)>` is invoked
with an initial step size such that the predicted plastic
strain delta does not exceed the limit.

:Type: real
:Default: :math:`-\infty`
:Valid Values: :math:`(-\infty,\infty)`

.. note::

   The legacy viscoplasticity solver defaulted this value to 1.1, which seems to
   moderately speed up calculations at the cost of stability.


nlk_tol (expert)
^^^^^^^^^^^^^^^^

The convergence tolerance for the NLK nonlinear solver for viscoplasticity. The
nonlinear system is considered solved by the current iterate if the norm of the
last solution correction is less than this value.

:Type: real
:Default: 1e-2
:Valid Values: (0, 1]


nlk_vector_tolerance (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vector drop tolerance for the NLK method. When assembling the acceleration
subspace vector by vector, a vector is dropped when the sine of the angle
between the vector and the subspace less than this value.

:Type: real
:Default: 0.01
:Valid Values: :math:`(0,1)`


pc_freq (expert)
^^^^^^^^^^^^^^^^

This controls how frequently the preconditioner is updated in the adaptive BDF2
integrator. A value of :math:`N` will allow a preconditioner to be used for as many
as :math:`N` consecutive time steps before being updated, although it may be updated
more frequently based on other criteria. A value of 1 causes the preconditioner
to be updated every time step.

:Type: integer
:Default: 1
:Valid Values: :math:`\geq 1`

.. note::

   A basic strategy of the adaptive BDF2 integrator is to use a preconditioner
   for as many time steps as possible, and only update it when a nonlinear time
   step iteration fails to converge. This generally works quite well. But if you
   find that the integrator is thrashing — evidenced by the number of times a
   step failed with an old preconditioner and was retried (this is the NNR
   diagnostic value in the terminal output) being a significant fraction of the
   number of time steps — it may be more cost effective to set this value to 1,
   for example.

solver (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The choice of viscoplastic solver. The default is fast and accurate, and an end
user won't benefit by changing this option. The default is to use the
``bdf2_integrator`` backend. The `"jacobian"` option uses an implicit
NLK-accelerated ``idaesol`` solver. The `"jfree"` option uses a jacobian-free
algorithm identical to the ``bdf2_integrator``, but implemented on the
``idaesol`` type.

:Type: string
:Default: `"bdf2"`
:Valid Values: `"bdf2"`, `"jacobian"`, or `"jfree"`

# Free-surface flow

`plug.json` is a reference-free, inviscid translating water plug. Water moves
at unit speed through a fluid/VOID pipe while the pressure boundary at the
upstream end is assigned `inflow-material: VOID`. Its regular 25×5 mesh leaves
full-fluid cells between the two interfaces, avoiding an under-resolved
two-interface reconstruction. The CTest regression checks the analytic plug
volume fractions, velocity, and pressure at three times. It is the
two-dimensional analogue of mainline `free-surf-flow-3`.

`plug_noisy.json` is the same plug problem on a seeded 10%-perturbed 25×5
mesh. Its analytic check clips the exact translated plug against the actual
cell polygons, testing nonorthogonal interface reconstruction without the
under-resolved two-interface configuration.

`fill.json` is the complementary water-inflow problem: water enters at unit
speed through a velocity boundary assigned `inflow-material: water`. It is the
two-dimensional analogue of mainline `free-surf-flow-1`. Together these cases
verify both VOID and material inflow routing to the geometric tracker.

`drain.json` ejects water through a unit-speed velocity outlet, while the
opposite pressure boundary is assigned `inflow-material: VOID`. It is the
two-dimensional analogue of mainline `free-surf-flow-2` and confirms the
matching analytic interface motion during drainage. The fill and drain cases
use the same 15×5 mesh, while the plug cases use 25×5.

`gravity_fill.json` and `gravity_drain.json` use zero pressure at both open
ends and an axial body acceleration instead of prescribed velocities. Their
continuum inviscid solutions have `u=t` and `u=-t`, respectively, with zero
pressure; the water/VOID interface follows `x=0.5+t²/2` for filling and
`x=4.5-t²/2` for draining. The tests use a small fixed step and compare the
interface with the corresponding exact trajectory of the implemented split
scheme, in which material transport uses the velocity at the beginning of the
step. They also check velocity, pressure, and water-volume conservation.

`hydrostatic.json` is a stationary water layer beneath a horizontal free
surface with downward body acceleration. It checks gravity-head balance,
zero pressure at the water/VOID interface, and the absence of spurious
velocity.

`hydrostatic_viscous.json` repeats the same case through the viscous momentum
path. Its physical solution should be identical, while the test also checks
that the viscous solver does not introduce spurious motion.

`hydrostatic_solid_wall.json` adds a pure SOLID strip along the left side of
the hydrostatic water/VOID configuration. It checks the interaction of
material-defined inactive cells, the free surface, and hydrostatic pressure
balance.

`input.json` is an exploratory inviscid fluid/VOID broken-dam problem in a
closed square domain. Water initially occupies the left half of the box and
evolves under downward gravity. It is the two-dimensional analogue of mainline
`free-surf-flow-4` and exercises geometric fluid/VOID tracking, free-surface
pressure treatment, and inactive-VOID flow cells.

# Free-surface flow

`plug.json` is a reference-free, inviscid translating water plug. Water moves
at unit speed through a fluid/VOID pipe while the pressure boundary at the
upstream end is assigned `inflow-material: VOID`. Its CTest regression checks
the analytic plug volume fractions, velocity, and pressure at three times. It
is the two-dimensional analogue of mainline `free-surf-flow-3`.

`fill.json` is the complementary water-inflow problem: water enters at unit
speed through a velocity boundary assigned `inflow-material: water`. It is the
two-dimensional analogue of mainline `free-surf-flow-1`. Together these cases
verify both VOID and material inflow routing to the geometric tracker.

`drain.json` ejects water through a unit-speed velocity outlet, while the
opposite pressure boundary is assigned `inflow-material: VOID`. It is the
two-dimensional analogue of mainline `free-surf-flow-2` and confirms the
matching analytic interface motion during drainage.

`input.json` is an exploratory inviscid fluid/VOID broken-dam problem in a
closed square domain. Water initially occupies the left half of the box and
evolves under downward gravity. It is the two-dimensional analogue of mainline
`free-surf-flow-4` and exercises geometric fluid/VOID tracking, free-surface
pressure treatment, and inactive-VOID flow cells.

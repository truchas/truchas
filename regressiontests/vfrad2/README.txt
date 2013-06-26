This test exercises the Skip_Geometry_Check=.true. feature of the
ENCLOSURE_RADIATION namelist.  The domain for this test consists of a
thin spherical shell containing a 'moving' ball.  Heat radiates between
the two bodies.  In the truchas simulation, HT is done with the ball
centered, but with view factors computed from a different mesh with the
ball shifted off center.  The truchas mesh and enclosure surfaces are
correctly linked in terms of faces, but the face geometries do not match
on the ball due to the shift.  The cell-based temperature results should
match those from the simulation that uses the off-center ball.

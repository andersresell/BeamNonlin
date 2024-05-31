# BeamNonlin
Project work for the PhD course KT8205.

A 3D explicit corotational beam element has been implemented. Since the formulation is explicit, only the inner forces are calculated, the tangent stiffness is left out. The implementation is based on the corotational beam in Crisfield: Nonlinear Finite Elements for Continua and Structures vol 2. The underlying locally linear beam element uses Euler Bernoulli beam thery to describe bending and standard axial and torsional rod thery to completely describe a 3D beam element. Nodal rotations are stored and updated using quaternions, and extracted as rotation matrices when needed in computations. Note that almost all published works on nonlinear beams uses implicit or quasistatic integration.

Since the resulting algorithm became unstable after some time for certain highly nonlinear configurations, a number of measures were attempted.
-   Both a rotation update using quaternions and the Hughes-Winget formula was tested, giving negligible differences.
-   Energy balance monitoring was implemented, confirming that the energy blew up in certain scenarios. The problem seems to be related to torsional vibrations. Reducing the time step     
    doesn't resolve the problem. 
-   A rotation time integration algorithm by Simo and Wong for rigid bodies has been used to integrate the rotations in time. This worked better than the algorithm for explicit integration 
    of rotations described in Crisfield's book (which was an extension of the half-step central difference for rotations). The latter method became unstable when trying to integrate a     
    freely spinning rigid body. The Simo and Wong method proved more robust, but didn't resolve the problem related to torsional vibrations.
-   The corotational Beam element of Battini was implemented in a debug attempt. The implementation was adopted from a python code found on github (not stated here for copyright reasons, 
    but google "corotational beam github" and it should appear). was verified that the inner forces produced were identical for the python code and for the code in this project (meaning         that my code has now been tested with inner force calculations from an existing (and hopefully verfified) code). The Battini beam formulation yielded very similar results as the 
    original Crisfield formulation, and the instability problems still persist. This measure strongly suggests that the issue is not related to the internal force calculations.

The reason for the instabilities are not yet discovered. Wether there is some bug lurking or if the corotational beam formulation is fundamentally problematic when used with explicit time integration is yet unknown to me.

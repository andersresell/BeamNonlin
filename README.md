# BeamNonlin
Project work for the PhD course KT8205.

A 3D explicit corotational beam element has been implemented. Since the formulation is explicit, only the inner forces are calculated, the tangent stiffness is left out. The implementation is based on the corotational beam in Crisfield: Nonlinear Finite Elements for Continua and Structures vol 2. The underlying locally linear beam element is described by Euler Bernoulli theory. Nodal rotations are stored and updated using quaternions, and extracted as rotation matrices when computing internal element forces. Time integration is handled with the momentum conserving Simo-Wong algorithm to solve rotations and the corresponding non-staggered central difference method to update translations at each node.

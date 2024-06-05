# BeamNonlin
Project work for the PhD course KT8205.

An explicit 3D beam solver based on a co-rotational formulation. Since the formulation is explicit, only the inner forces are calculated, the tangent stiffness is left out. The implementation is based on the corotational beam by Crisfield, the formulation can be found in "Nonlinear Finite Elements for Continua and Structures vol 2". The underlying locally linear beam element is described by Euler Bernoulli theory. Nodal rotations are stored and updated using quaternions, and extracted as rotation matrices when computing internal element forces. Time integration is handled with the momentum conserving Simo-Wong algorithm to solve rotations and the corresponding non-staggered central difference method to update translations at each node.
![example](https://raw.githubusercontent.com/andersresell/BeamNonlin/main/beam-showcase.gif)

Dependencies:
- Linux enviroment
- build-essential
- eigen
- yaml-cpp
- To plot: python with numpy and matplotlib

To test the code, clone it and add: export BeamNonlinHome=/path/to/source/directory to ~/.bashrc
The project is built and a default test is ran by running ./run-test.sh from the source directory.
The code can be used from a different directory by running $BeamNonlinHome/launch.sh <input-file.yml> where input-file.yml is an input file specifying the problem.

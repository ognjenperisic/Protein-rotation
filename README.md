# Protein-rotation

The complexity of classical molecular dynamics simulations linearly depends on the number of simulated atoms (doubling the number of atoms doubles the CPU load/simulation time). 


The first code depicted here, **rotation.cpp**, reorients a given molecule (pdb file), so its longest dimension is aligned with one of the principal axes (X, Y or Z). That reduces the number of atoms in a simulation box (hint: Pythagorean Theorem). 


The second code, **rotation_2.cpp**, reorients a molecule (pdb file), so the vector between two given residues is aligned with the X axis. This is useful in steered molecular dynamics simulations. The similar functionality exists in most major molecular modelling tools.


The code is written in bare bones C++, with no additional libraries. 

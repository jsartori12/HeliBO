# Helical Binder Optimizer (HeliBO)
# Description
HeliBO is a python script that uses the pyRosetta software to optimize the binding free energy of helical peptides to a specific target. This is done by iteratively optimizing the peptide using the Rosetta software's Design function. At each optimization round, the peptide is selected according to the criteria defined by the user, who can choose between local search and Monte Carlo. In addition, the script has the option of activating a composition constraint in order to prioritize the selection of residues with a high propensity to form an α-helix using AddHelixSequenceConstraints Mover, which can be used after applying a weight to the aa_composition energy term. The user can define the weight they want to give to the aa_composition term, where higher weights try to further prioritize the selection of those residues with a high propensity to form an α-helix.


# Usage

python Optimization_cycle.py --pdb  --chain  --unbind  --constraints  --weight  --decision  --n_cycle 

--pdb = PDB structure of the protein to optimize

--chain = Chain ID to optimize

--unbind = Chains selecteds to unbind and calculates de dG values
  - usage: A_D, to unbind chains A and D 

--constraints = True or False
  - If True = Add aa_composition term in energy function and apply HelixConstraintMover.
  - If False = Don't use aa_composition term
  
--weight = 0 - 1.0
  - If constraints = True, apply the input weight aa_composition term.

--decision = MC or LS
  - If MC = calculates ddG and apply Metropolis Monte Carlo criteria make the decision
  - If LS = calculates ddG and apply Local Search in the decision

--n_cycle = number of rounds for the optimization cycle

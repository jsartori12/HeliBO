# Peptide-Optimization-Cycle



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

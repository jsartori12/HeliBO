# Helical Binder Optimizer (HeliBO)
# Description
HeliBO is a python script that uses the pyRosetta software to optimize the binding free energy of helical peptides to a specific target. This is done by iteratively optimizing the peptide using the Rosetta software's Design function. At each optimization round, the peptide is selected according to the criteria defined by the user, who can choose between local search and Monte Carlo. In addition, the script has the option of activating a composition constraint in order to prioritize the selection of residues with a high propensity to form an α-helix using AddHelixSequenceConstraints Mover, which can be used after applying a weight to the aa_composition energy term. The user can define the weight they want to give to the aa_composition term, where higher weights try to further prioritize the selection of those residues with a high propensity to form an α-helix.




# Installing dependencies

## Getting started

First of all, you must download PyRosetta. To download , you need to get a license.
<br />
License and Downloads links:
<br />
[License](https://www.rosettacommons.org/software/license-and-download)
<br />
[PyRosetta4 Download](https://graylab.jhu.edu/download/PyRosetta4/archive/release/)



## Installing PyRosetta
### After downloading, unzip PyRosetta's and enter the setup directory to install it
```
tar -xvf PyRosetta[release].tar.bz2
cd PyRosetta*/setup
python3 setup.py install
```
## Installing aditional Python libs
```
pip install pandas
pip install numpy
```

Additional help for downloading and installing and PyRosetta (source:Sari Sabban youtube channel )

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/UEaFmUMEL9c/0.jpg)](https://www.youtube.com/watch?v=UEaFmUMEL9c)


## Download the repository in your directory


```
cd your_directory
git clone https://github.com/jsartori12/HeliBO.git
```



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

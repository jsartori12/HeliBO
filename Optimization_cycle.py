#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 20:00:41 2023

@author: joao
"""

import argparse
import pandas as pd
import os
import random
import math
#imports from pyrosetta
from mimetypes import init
from pyrosetta import *
from pyrosetta.teaching import *
#from IPython.display import Image
#Core Includes
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols import *

pyrosetta.init()
#pyrosetta.init(options="-constant_seed -jran 314")

scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

#### create Dumps directory to save pdbs and CSV

#new_dir = "Dumps"
#ac_dir = "./"

#path = os.path.join(ac_dir, new_dir)
#os.mkdir(path)

#path_new = os.path.join(path, "PDBs")
#os.mkdir(path_new)

####
#### FastRelax Protocol - bb_unlocked
####
def pack_relax(pose, scorefxn):
    """Function that performs FastRelax protocol from pyRosetta.
    Parameters
        ----------
        pose : pose_obj
            pose object from pyrosetta
        scorefxn : score_obj
            score function defined for a given pose"""
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    # Set up a MoveMapFactory
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)

    ## Print informations about structure before apply fast relax
    # display_pose = pyrosetta.rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
    # display_pose.tasks(tf)
    # display_pose.movemap_factory(mmf)
    # display_pose.apply(pose)

    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    return 

#### Apply design to inputed chain in a Pose - 
####
def design(pose, chain, constraint):
    """Function used to design protein. Use design function from pyRosetta
    to select the best romater in each position for the selected chain. 
    Parameters
        ----------
        pose : pose_obj
            pose from protein to design
        chain : str
            chain ID from the protein to apply the design function
        constraint : bool
            compositional constraint mover to apply AddHelixSequenceConstraintsMover.
            Mover favors helix-prone residues during the design.
            """
    constraint = args.constraints
    if constraint == False:
        #### Get pose residues for selected chain to allow design in these positions
        chain_D = pyrosetta.rosetta.core.pose.get_resnums_for_chain(pose, chain)
        chain_D = list(chain_D)

        mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        mut_posi.set_index(str(chain_D[0])+"-"+str(chain_D[-1]))
        #### Print selected residues for mut posi
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))

        # Select Neighbor Position
        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(mut_posi)
        nbr_selector.set_include_focus_in_subset(True)
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

        # Select No Design Area
        not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
        #### Print residues to NOT design
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

        # The task factory accepts all the task operations
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

        # These are pretty standard
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

        # Disable Packing
        prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
        tf.push_back(prevent_subset_repacking)

        # Disable design
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

        # Enable design
        aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        aa_to_design.aas_to_keep("ACDEFHIKLMNQRSTVWY")
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

        # Create Packer
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
        packer.task_factory(tf)
        # Fast relax protocol
        mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
        mmf.all_bb(setting=True)
        mmf.all_bondangles(setting=True)
        mmf.all_bondlengths(setting=True)
        mmf.all_chi(setting=True)
        mmf.all_jumps(setting=True)
        mmf.set_cartesian(setting=True)
        # Apply fastrelax to pose
        fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
        fr.cartesian(True)
        fr.set_task_factory(tf)
        fr.set_movemap_factory(mmf)
        fr.min_type("lbfgs_armijo_nonmonotone")
        ## Apply design to pose
        packer.apply(pose)
        fr.apply(pose)
        return pose
    #### If Constraints set to True, apply a weigth to aa_composition energy term and also apply AddHelixSequenceConstraintsMove
    #### to the design function.
    else:
        #### Get pose residues for selected chain to allow design in these positions
        chain_D = pyrosetta.rosetta.core.pose.get_resnums_for_chain(pose, chain)
        chain_D = list(chain_D)

        mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        mut_posi.set_index(str(chain_D[0])+"-"+str(chain_D[-1]))
        #### Print selected residues for mut posi
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))

        # Select Neighbor Position
        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(mut_posi)
        nbr_selector.set_include_focus_in_subset(True)
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

        # Select No Design Area
        not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
        #### Print residues to NOT design
        #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

        # The task factory accepts all the task operations
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

        # These are pretty standard
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

        # Disable Packing
        prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
        tf.push_back(prevent_subset_repacking)

        # Disable design
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

        # Enable design
        aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        aa_to_design.aas_to_keep("ACDEFHIKLMNQRSTVWY")
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

        # # Clear constraints
        clear = pyrosetta.rosetta.protocols.aa_composition.ClearCompositionConstraintsMover()
        clear.apply(pose)
        #### Add helix constraints and assign a weight to the aa_composition term (weight defined in the input script)
        #### Information about the mover in: 
        #### https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/best_practices/AddHelixSequenceConstraintsMover
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.aa_composition, args.weight)
        add_helix_sequence_constraints = pyrosetta.rosetta.protocols.aa_composition.AddHelixSequenceConstraintsMover()
        add_helix_sequence_constraints.set_residue_selector(mut_posi)
        add_helix_sequence_constraints.apply(pose)
        # print(scorefxn.show(pose))
        # Create Packer
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
        packer.task_factory(tf)
        # Fast relax protocol
        mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
        mmf.all_bb(setting=True)
        mmf.all_bondangles(setting=True)
        mmf.all_bondlengths(setting=True)
        mmf.all_chi(setting=True)
        mmf.all_jumps(setting=True)
        mmf.set_cartesian(setting=True)
        # Apply fastrelax to pose
        fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
        fr.cartesian(True)
        fr.set_task_factory(tf)
        fr.set_movemap_factory(mmf)
        fr.min_type("lbfgs_armijo_nonmonotone")
        ## Apply design to pose
        packer.apply(pose)
        fr.apply(pose)
        return pose
#### Unbind chains in pose
####
def unbind(pose, partners):
    """Function that performs unbind of complex to allow dG calculation.
    Parameters
        ----------
        pose : pose_obj
            pose object from pyrosetta
        partners : str
            name of chains to unbind. usage example to split chains A and B from a protein complex: 'A_B'  """
    #### Generates dummy pose to maintain original pose
    pose_dummy = pose.clone()
    pose_binded = pose.clone()
    STEP_SIZE = 100
    JUMP = 1
    docking.setup_foldtree(pose_dummy, partners, Vector1([-1,-1,-1]))
    trans_mover = rigid.RigidBodyTransMover(pose_dummy,JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose_dummy)
    pack_relax(pose_dummy, scorefxn)
    #### Return a tuple containing:
    #### Pose binded = [0] | Pose separated = [1]
    return pose_binded , pose_dummy

#### Calculates dGbinding - Unbounded - bounded dG
####
def dG_v2_0(pose_Sep, pose_bind):
    """Function that performs dG calculation for a protein complex.
    Parameters
        ----------
        pose_Sep : pose_obj
            pose object with chains separated after applying unbind function
        pose_bind : pose_obj
            pose object with chains bounded  """
    bound_score = scorefxn(pose_bind)
    unbound_score = scorefxn(pose_Sep)
    dG = bound_score - unbound_score
    return dG

#### Decision, calculates ddG - Monte Carlo or Local Search
####
def decision(before_pose, after_pose, dG_before, dG_after, decision_method):
    """Function used to apply decision criteria. Compare a structure before and
    after design and decided which based on dG value.
    Parameters
        ----------
        before_pose : pose_obj
            pose from protein before apply design function
        after_pose : pose_obj
            pose from protein after apply design function 
        dG_before : double
            dG value from protein before apply design function
        dG_after : double
            dG value from protein after apply design function
        decision_method : str
            decides which decision criteria to apply.
            Criterias available: Local Search, Monte Carlo.
            Usage: 'MC' for Monte Carlo; 'LS' for Local Search
            """
    #### Calculates ddG and pick    T H E   C H O S E N   O N E
    decision_method = args.decision
    E = dG_after - dG_before
    if decision_method == "MC":
        if E < 0:
            return after_pose
        elif random.uniform(0, 1) >= math.exp(-E/1):
            return before_pose
        else:
            return after_pose
    elif decision_method == "LS":
        if E < 0:
            return after_pose
        else:
            return before_pose

#### Optimization cycle
####   
def Optimization_constraint(pose,n_cycle, constraints):
    """Function used to apply decision criteria. Compare a structure before and
    after design and decided which based on dG value.
    Parameters
        ----------
        pose : pose_obj
            pose from protein to design
        n_cycle : int
            number of cycle for optimization
        constraint : bool
            compositional constraint mover to apply AddHelixSequenceConstraintsMover.
            Mover favors helix-prone residues during the design.
            """
    n_cycle = args.n_cycle
    constraints = args.constraints
    deltasGes = []
    peptideenergy = []
    sequences = []
    if constraints == True:
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.aa_composition, args.weight)
    pose_before = unbind(pose, args.unbind)
    sequences.append(pose_before[1].sequence())
    deltasGes.append(dG_v2_0(pose_before[1], pose_before[0]))
    for i in range(n_cycle):
        pose_before = pose_before
        pose_before_temp = pose_before[0].clone()
        #### Apply Design
        after_pose = design(pose_before_temp, args.chain, args.constraints)
        #### Unbind designed pose
        after_pose = unbind(after_pose, args.unbind)
        #### Compares pose before and after design and returns THE CHOSEN ONE
        decided_pose = decision(pose_before, after_pose,dG_v2_0(pose_before[1], pose_before[0]),dG_v2_0(after_pose[1], after_pose[0]), args.decision)
        #### Appends THE CHOSEN ONE sequence and dG
        sequences.append(decided_pose[0].sequence())
        deltasGes.append(dG_v2_0(decided_pose[1], decided_pose[0]))
        #### Save THE CHOSEN ONE pose as a PDB
        decided_pose[0].dump_pdb("./Dumps/PDBs/log_mutated_"+str(i+1)+".pdb")
        #### Save data in a .csv file
        data = {"dG: ": deltasGes,
            "Sequences": sequences}
        df = pd.DataFrame(data)
        df.to_csv("./Dumps/teste.csv")
        #### Inserts THE CHOSEN ONE pose into pose_before to return de cycle
        pose_before = decided_pose
if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--chain', type=str, required=True)
    parser.add_argument('--unbind', type=str, required=True)
    parser.add_argument('--constraints', type=bool, required=True)
    parser.add_argument('--weight', type=float, required=False)
    parser.add_argument('--decision', type=str, required=True)
    parser.add_argument('--n_cycle', type=int, required=True)


    args = parser.parse_args()


    pose = pose_from_pdb(args.pdb)
    scorefxn.score(pose)
    Optimization_constraint(pose, args.n_cycle, args.constraints)

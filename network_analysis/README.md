# Network analysis

This module is meant for network analysis of created models. That means that the network properties of the model are analyzed without biological information.

The submodules include:
- find_components.py
- loop_removal.py
- compaction.py

## Find Components

The find components submodule identifies and visualizes components of the module network that are unconnected to the rest of the network.

## Loop Removal

The loop removal module identifies reactions that are only involved in thermodynamically infeasibile loops. These reactions are obligatory loop reactions meaning they can only participate in infeasible loops. It does not identify reactions that can conditionally participate in loops. Therefore, it means that the module doesn't remove all of the possible loops from the model. 

The loop removal uses the compaction algorithm to identify obligatory loop reactions.

## Compaction

The compaction algorithm compacts linear pathways and combines identical reactions. The linear compaction is done by identifying metabolites that participate only in 2 reactions. Such metabolites can be removed by combining the 2 linked reactions appropriately. Combining identical reactions is done by comparing reactions. If 2 reactions have the same stoichiometric equation or can be scaled to have the same equation, they are combined. Both of the linear compaction and reaction combining is done for several iterations until no more compaction actions can be identified.

Compacted reactions that have no metabolites after compaction correspond to thermodynamically infeasible loops. In such such cases, the compacted reactions cancel out and show that the Gibbs free energy is 0 for such cycle. The simplest example is 2 identical reactions that have the opposite directions:

        ┌──── (R1) ────►
  ---> [A]             [B]
        └──── (R2) ◄────┘

Althogh R1 and R2 are connected to the rest of the network through A, they are still obligatory loop reactions as B can only be produced and consumed in a loop. If R1 and R2 are replaced by linear or parallel networks, a less trivial example can be imagined.
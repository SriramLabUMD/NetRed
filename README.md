# NetRed
NetRed, an algorithm to reduce genome-scale metabolic networks and facilitate the analysis of flux predictions

## Publication

**If used, please cite the following article**

Lugar, D., Mack, S., & Sriram, G. (2020). NetRed, an algorithm to reduce genome-scale metabolic networks and facilitate the analysis of flux predictions. *Metabolic Engineering*, doi: 10.1016/j.ymben.2020.11.003


## Instructions

This program is intended for use in reducing the results of Flux Balance Analysis (FBA) simulations. For details of the algorithm, please refer to the associated publication.

### Inputs 

**NetRed(S,fluxes,mets,protectedMets,options)**  

#### S
The stoichiometric matrix  
  
#### fluxes
1. If analyzing a single flux solution: Enter a single flux (column) vector solution corresponding to the given stoichiometric matrix, with the fluxes in the same order as their corresponding columns in S.  
2. If comparing multiple flux solutions: Enter a matrix whose columns are composed of the flux solution vectors being analyzed.  
               
#### mets
A complete list of metabolite names in the same order as their corresponding rows in S.  

#### protectedMets
A list of metabolite names desired to remain in the network by the user after network reduction.  

#### options
User specified options may be included using the following fields:  

*'filename'* - Name of the file to which the results are written  
*'threshold'* - Minimum value at which a flux is considered active (default 1E-8)  
*'lumpIdenticalRxns'* - Flag (1 or 0) to combine parallel reactions into one after network reduction (default 0)  
*'verbose'* - Flag (1 or 0) to print algorithmic details while reducing the network (default 1)  
*'BM'* - Flag (1 or 0) to resolve degenerate biomass equations  
         

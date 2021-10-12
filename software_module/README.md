## `wgregseq`

This is a python software module which includes all custom written code needed for this work. The module can be installed locally after downloading or cloning the repository. Therefore navigate into this folder in the command line on your machine and type

`pip install -e .`

Afterwards the package can simply be imported in any python code. 


### Modules
There are various submodules that contain functions similar contexts. 

 - `viz.py` : Code used to visualize results at various stages of the pipeline.
 - `seq_utils.py` : Functions that create and/or manipulate sequences. 
 - `model_utils.py` : Any helper functions to use thermodynamic modeling, e.g., creating and evaluating energy matrices.
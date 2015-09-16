# alignedSENSE
Tools for aligned reconstruction of multishot MR.

This repository provides tools to implement the reconstruction methods and reproduce the synthetic experiments included in the manuscript ''Sensitivity encoding for aligned multishot magnetic resonance reconstruction'', L. Cordero-Grande, R. P. A. G. Teixeira, E. J. Hughes, J. Hutter, A. N. Price, and J. V. Hajnal. Unpublished.

The code has been developed in MATLAB and has the following structure:

###### ./
contains the scripts for running the experiments included in Section III of the manuscript: *alignedSENSE_Exp[1-4].m* and generate the corresponding convergence plots: *plot_Exp[1-4].m*.

###### ./data
contains the dataset used for simulations: *xGT.mat*. Data generated when running the scripts is also stored in this folder.

###### ./lib
contains a package used to generate plot colouring: *varycolor*.

###### ./meth
contains the solvers for motion and reconstruction described in the manuscript: *solve[X,T].m*.

###### ./meth/pre
contains functions used for preprocessing: *coilArrayCompression.m*.

###### ./meth/pro
contains functions used by the solvers: *fftGPU.m*, *ifftGPU.m*, *sense.m*, *isense.m*, *precomputationsSincRigidTransform.m*, *sincRigidTransform.m*, *sincRigidTransform[Gradient,Hessian].m*.

###### ./meth/pos
contains the main function used to evaluate the results: *errorFit.m*.

###### ./synth
contains the methods used to generate synthetically motion corrupted data: *generate[Encoding,Grids].m*, *synthesize[Y-T].m*.

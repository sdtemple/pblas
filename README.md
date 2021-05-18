# pblas
My implementations for pair-based likelihood approximations as seen in Stockdale, et al. (2019). This software adds to the original reference code at jessicastockdale/PBLA. Namely, I have:

* implemented all PBLAs in the R scripting language;
* optimized a O(m^2 n^2) task for computational efficiency
* provided wrappers `pbla_gsem` and `pbla_multi` for optimizing general and multitype stochastic epidemic models;
* made the option of a fixed exposed period for SEIR models;
* made the option of numerous patient zeros;
* and supplemented code for simulation studies, MCMC samplers, and real data analyses.

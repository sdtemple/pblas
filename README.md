# pblas
My implementations for pair-based likelihood approximations as seen in Stockdale, et al. (2019). This software adds to the original reference code at [jessicastockdale/PBLA](https://github.com/jessicastockdale/PBLA). Namely, I have:

* implemented all PBLAs in the R scripting language;
* optimized an O(m^2 n^2) task for computational efficiency;
* provided wrappers `pbla_gsem` and `pbla_multi` for optimizing general and multitype stochastic epidemic models (SEMs);
* made the option of a fixed exposed period for SEIR models;
* made the option of numerous patient zeros;
* and supplemented code for simulation studies, MCMC samplers, and real data analyses.

Some functions exist to reproduct the results of Stockdale, et al. (2019). Practitioners should use:

* `pbla_std` for complex SEMs;
* `pbla_prod` for general SEMs;
* `rgsem` to simulate epidemics;

For large epidemics, some efficient memory use may be important; I demonstrate such with global variable assignments `<<-` in `pbla-ebola.R`. To simulate a general stochastic epidemic, use `rgsem`. (See `multitypes.R` with `pbla_multi.R` for mulitype SEM. A case study is `pbla-mcmc-tristan.R`.)

Some key files in this repository are:

* `pbla-report-sdtemple.pdf` documents the model, theoretical developments, and experimental results.
* `pbla-slides-sdtemple.pdf` summarizes the model, theoretical developments, and experimental results.
* `pbla-bakeoff.R` compares 5 PBLAs in a simulation study.
* `pbla-ebola.R` examines the Ebola virus epidemic in West Africa.
* `pbla-nN.R` studies PBLA independence assumptions under increasing n / N infected proportion. 
* `pbla-mcmc-tristan.R` examines the common cold outbreak in Tristan da Cunha.
* `pbla-exp.R` is the basis for simulation studies with exponential infectious periods.
* `pbla-erl.R` is the basis for simulation studies with Erlang infectious periods.
* `pbla-mcmc-example.R` offers a sampler with random walk proposals and Hastings ratios based on PBLAs.
* `ncda-mcmc-example.R` offers a sampler for non-centered data augmentation.
* A sample script for simulation study figures is `pbla-figure.R`. I modify this script to generate various figures. 

This package has been tested on R 3.6.2 and R 4.0.5. Besides `pbla_prod_parallel`, source code is exclusively base R, so the package should operate for older versions of R. `pbla_prod_parallel` utilizes `parallel`, `foreach`, and `doParallel` to offer parallel computing for the O(n^2) loop. Parallel computing will only speed up calculations for n > 10,000. 

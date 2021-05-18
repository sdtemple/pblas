# pblas
My implementations for pair-based likelihood approximations as seen in Stockdale, et al. (2019). This software adds to the original reference code at [jessicastockdale/PBLA](https://github.com/jessicastockdale/PBLA). Namely, I have:

* implemented all PBLAs in the R scripting language;
* optimized a O(m^2 n^2) task for computational efficiency
* provided wrappers `pbla_gsem` and `pbla_multi` for optimizing general and multitype stochastic epidemic models;
* made the option of a fixed exposed period for SEIR models;
* made the option of numerous patient zeros;
* and supplemented code for simulation studies, MCMC samplers, and real data analyses.

Some key files in this repository are:

* `pbla-report-sdtemple.pdf` documents the model, theoretical developments, and experimental results.
* `pbla-slides-sdtemple.pdf` summarizes the model, theoretical developments, and experimental results.
* `pbla-bakeoff.R` compares 5 PBLAs in a simulation study.
* `pbla-ebola.R` examines the Ebola virus epidemic in West Africa.
* `pbla-mcmc-tristan.R` examines the common cold outbreak in Tristan da Cunha.
* `pbla-exp.R` is the basis for simulation studies with exponential infectious periods.
* `pbla-erl.R` is the basis for simulation studies with Erlang infectious periods.
* `pbla-mcmc-example.R` offers an example MCMC sampler with random walk proposals and Hastings ratios based on PBLAs.
* `ncda-mcmc-example.R` offers and example MCMC sampler for non-centered data augmentation.
* Various sample scripts for figures in my report are under `figures`. 

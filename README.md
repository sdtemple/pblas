# pblas
My implementations for pair-based likelihood approximations as seen in Stockdale, et al. (2019). This software adds to the original reference code at [jessicastockdale/PBLA](https://github.com/jessicastockdale/PBLA). Namely, I have:

* implemented all PBLAs in the R scripting language;
* optimized an O(m^2 n^2) task for computational efficiency;
* provided wrappers `pbla_gsem` and `pbla_multi` for optimizing general and multitype stochastic epidemic models;
* made the option of a fixed exposed period for SEIR models;
* made the option of numerous patient zeros;
* and supplemented code for simulation studies, MCMC samplers, and real data analyses.

See [this video](https://youtu.be/jv1vYEU-VNA) for my oral exam. Some functions exist to reproduct the results of Stockdale, et al. (2019). Practitioners should use:

* `pbla_std` for complex SEMs;
* `pbla_prod` for general SEMs;
* `rgsem` to simulate epidemics;
* `pbla_gsem`, `pbla_multi` as wrappers for MLE.

For large epidemics, some efficient memory use may be important; I demonstrate such with global variable assignments `<<-` in `pbla-ebola.R`. To simulate a general stochastic epidemic, use `rgsem`. (See `multitypes.R` with `pbla_multi.R` for mulitype SEM. A case study is `pbla-mcmc-tristan.R`.)

Some key files in this repository are:

* `pbla-report-sdtemple.pdf` documents the model, theoretical developments, and experimental results;
* `pbla-slides-sdtemple.pdf` summarizes the model, theoretical developments, and experimental results;
* `pbla-ebola.R`, `pbla-ebola-contours.R` examine the Ebola virus epidemic in West Africa (see `data/ebola/`);
* `pbla-nN.R` studies PBLA independence assumptions under increasing n / N infected proportion;
* `pbla-under.R` considers adjustments in the case of underreporting;
* `pbla-time.R`, `pbla-parallel.R` benchmark computational performance;
* `pbla-rabies` investigates the dog rabies epidemic in Bangui, Central African Republic; 
* `pbla-mcmc-tristan.R` examines the common cold outbreak in Tristan da Cunha;
* `pbla-mcmc-example.R` offers a sampler with random walk proposals and Hastings ratios based on PBLAs;
* `ncda-mcmc-example.R` offers a sampler for non-centered data augmentation;
* `pbla-bakeoff.R` compares 5 PBLAs in a simulation study;
* `pbla-exp.R` is the basis for simulation studies with exponential infectious periods;
* `pbla-erl.R` is the basis for simulation studies with Erlang infectious periods;
* `pbla-figure.R` is a sample script for figures in simulation studies. 

This package has been tested on R 3.6.2 and R 4.0.5. Besides `pbla_prod_parallel`, source code is exclusively base R, so the package should operate for older versions of R. `pbla_prod_parallel` utilizes `parallel`, `foreach`, and `doParallel` to offer parallel computing for the O(n^2) loop. Parallel computing will only speed up calculations for n > 10,000.

Correction! Please not that there is a mistake in the slides and in the report. On slide 8/31, the MLE for the removal rate gamma should be the MLE for an exponential rate (lambda on Wikipedia article for exponential distribution). This impacts the figure on slide 22/31. R[0] estimates for the completely observed SEM do increase with increasing infected proportion. Thus, any commentary in the report and in the YouTube talk are invalid.

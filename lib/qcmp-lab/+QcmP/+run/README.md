# `QcmP.run`
Running [QcmP-lab](https://github.com/QcmPlab) codes from MATLAB

----------

### Choose your favorite workflow


#### Production

- `dry_line()` The most basic scenario: linear increase of U, no feedback, just some noninvasive check of convergence to flag out, _a posteriori_, which points should be discarded. **[deprecated]**

- `vec_line()` New API for dry lines (see above): instead of passing a start, step and stop values for a linear span, you just directly pass an array, built the way you prefer. The original behavior of `dry_line()` is obtained by simply passing `start:step:stop`, as easy as it gets. Silly me for not writing it this way the firs time.

- `autostop_line()` Basic linspace in U, without any adaptive feedback mechanism. Just stops the moment dmft does not converge (to avoid wasting cpu-time).

- `autostep_line()` Gradually increases U by feedback-controlled steps, i.e. if dmft does not converge the point is discarded and the step reduced.


#### Refinement

- `refresh_line()` Systematically enters the folder structure of a pre-existent calculation and uses the given _used_ files as a _restart_ to perform some additional dmft loops. Useful if something has been added to the driver or if you want to change something in the inputfile for just one last loop (e.g. `dm_flag`, `chi_flag`...). It has also a single-point variant, for maximum flexibility: `refresh_point()` that can be called from within a computation folder.

- `inspect_line()` Call if you wish to inspect, at runtime, how much points have been computed, how much have converged and which have not.


#### Testing & Setup

- `single_point()` It runs just a single-point calculation, but it supports a `Uold` value, to handle restart-files; to be called from a blank folder (the `inputfile` should be there, `Uold=%f` folder too).

----------

### TODO

- [ ] `dry_array()` A variant of `dry_line()` that exploits the array-env-variables provided by `SLURM` on HPC facilities.

- [ ] `refresh_array()` A variant of `refresh_line()` that exploits the array-env-variables provided by `SLURM` on HPC facilities.

- [ ] `autostep_array()` A variant of `autostep_line()` that exploits the array-env-variables provided by `SLURM` on HPC facilities.

- [ ] `autosampling_line` A replacement to `autostep_line()` (which would enter a deprecation cycle), aimed at controlling the step by inspecting the evolution of the number of loops.

### NOT-TODO

- [ ] `automixing_line()` Automatic control of self-mixing. Requires inspecting evolution of the dmft error, could be cumbersome (better handled at a lower level, e.g. within the dmft-loop driver: see [`adaptive_mix()`](https://github.com/QcmPlab/SciFortran/blob/master/src/SF_OPTIMIZE/adaptive_mix.f90) in SciFortran).

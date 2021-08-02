# treeducken 1.1.1

## Major changes

* Added mutualism boolean for simulation
* Switch `hs_mode` from a boolean to a string. Can now use `switch`, `spread`, or `both`.


## Bug fixes

* Fixed a bug where instead of skipping a replicate when using `summarize_cophy` the function would stop when encountering a dataset where the Parafit couldn't be calculated.
* Fixed a bug where `summarize_1cophy` would have NA as the column name if the column value had NA


# treeducken 1.1.0

## Major changes

* Function names shortened to be more user friendly
* Added functionality for plotting events on cophylogeny from `event_history` data frame
* `sim_cophyloBD` (formerly `sim_cophylo_bdp`) now has two arguments:
    * `host_limit` controlling the number of hosts a symbiont is allowed to have at one time
    * `hs_mode` which changes host expansion events to host switch events (`h_exp_r` is still used to change the rate)
* Added new simulation mode `sim_cophyloBD_ana` which adds anagenetic events (i.e., happening between speciation nodes)  to the functionality of `sim_cophyloBD`
* Changed event history to output only major events (description of abbreviations added to documentation)
* `show_div_bar` option added to plot. This outputs a bar showing the different events as tick marks.
* swapped output of `association_mat` from hosts in columns and symbionts in rows to hosts in rows and symbionts in columns.

# Bug fixes

* Fixed typo in README.md


# treeducken 1.0.0

## Major changes

* First release of software can currently simulate:
    * species trees
    * host-symbiont paired trees
    * locus trees within host, symbiont, or species trees
    * gene trees within locus trees using the multilocus coalescent
    * gene trees within species trees using the multispecies coalescent

* Other features:
    * plotting of host-symbiont paired trees that include extinct tips
    * limited functions for using the simulated trees

## Bug fixes

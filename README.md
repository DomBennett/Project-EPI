# Quantifying the living fossil with an Evolutionary Performance Index (EPI)

This pipeline generates EPI values for any group provided a tree and character states are
provided.

*Scripts*

* `run.R`: run all stages with given parameters
* `parameters.R`: set parameters for pipeline (open script, set parameters, run script)
* `setup.R`: ensure data folder is in order, ensure dependencies are installed

*Stages*

0. `0_data.R`: wrangle downloaded data from 0_data/ for subsequent stages
1. `1_change.R`: estimate change
2. `2_time.R`: estimate time
3. `3_success.R`: estimate success
4. `4_epi.R`: calculate EPI
5. `5_correlates.R`: determine correlates

*Reference*

**Not yet published**

*Author*

Dom Bennett
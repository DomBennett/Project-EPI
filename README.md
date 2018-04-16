# Quantifying the living fossil with an Evolutionary Performance Index (EPI)

Towards a quantified living fossil definition through the development of an evolutionary performance index (EPI). This pipeline generates EPI values for (potentially) all clades across the metazoa and plants.
It requires the following initial data files to be downloaded and placed in the `0_data/` folder:
[NCBI taxonomy .dmp](ftp://ftp.ncbi.nih.gov/pub/taxonomy/), [panTheria database](http://esapubs.org/archive/ecol/e090/184/), [O'Leary mammalian morpholopies](http://www.morphobank.org/index.php/Projects/ProjectOverview/project_id/773),
[Livezy avian morphologies](http://www.ncbi.nlm.nih.gov/pubmed/18784798) and [Lislevand avian traits](http://www.esapubs.org/archive/ecol/E088/096/default.htm#data). (URLs accessed Sept. 2016).

**Pipeline Structure**

![EPI pipeline](other/pipeline.png?raw=true "EPI pipeline")


**Folder Organisation**

```
-- Project-EPI/
-- -- run.R
-- -- parameters.R
-- -- setup.R
-- -- 0_data/
-- -- -- chars
-- -- -- ncbi_taxonomy
-- -- -- raw
-- -- -- trees
-- -- 1_wrngl/
-- -- 2_chng/
-- -- 3_node_obj/
-- -- 4_cntrst_n/
-- -- 5_split/
-- -- 6_phylotime/
-- -- 7_timetree/
-- -- 8_cntrst_chng/
-- -- 9_epi/
-- -- stages/
-- -- -- 1_wrngl.R
-- -- -- 2_chng.R
-- -- -- 3_node_obj.R
-- -- -- 4_cntrst_n.R
-- -- -- 5_split.R
-- -- -- 6_phylotime.R
-- -- -- 7_timetree.R
-- -- -- 8_cntrst_chng.R
-- -- -- 9_epi.R
-- -- -- analysis.R
-- -- tools/
-- -- -- analysis_tools.R
-- -- -- chng_tools.R
-- -- -- clade_matching_tools.R
-- -- -- epi_tools.R
-- -- -- i_tools.R
-- -- -- node_obj_tools.R
-- -- -- phyloprmttn_tools.R
-- -- -- phylotime_tools.R
-- -- -- timetree_tools.R
-- -- -- wrngl_tools.R
-- -- caches/
-- -- other/
```

`run.R` will run all the stages by loading the parameters from `parameters.R`. Edit `parameters.R` to change the way the pipeline works. All stage scripts are held in a separate folder, the specific functions they requier to run are held in `tools/`. For a quick-and-dirt setup use `setup.R` to automatically install dependencies and ensure the folders are in the right structure with correct initial files (this script is not infallible). `other/` contains leftover scripts and explorations during pipeline development, `caches/` contains data downloaded from [timetree](http://www.timetree.org/), by storing it locally fewer request are required.


**Running**

The entire pipeline can be run via the setup and run scripts, e.g. in UNIX:

```{bash}
Rscript setup.R
Rscript run.R >& log &
```

**Reference**

Bennett, Dominic J., Sutton, Mark D., and Turvey, Samuel T. 2018. Quantifying the living fossil concept. *Palaeontologia Electronica* 21.1.15A 1-25. https://doi.org/10.26879/750

**Author(s)**

[Dom Bennett](https://github.com/DomBennett)

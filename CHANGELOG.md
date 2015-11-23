# TRENTOOL CHANGE LOG

## Unreleased
- new plotting features in ```TEplot2D_beta.m```

## 3.4.2 - 2015-11-23

### Fixed
- removed unnecessary output variables, closes [#16]
- bug fix in ```TEgroup_stats.m```, unordered data are now collected correctly for dependent samples tests

### Added
- Statistics on TE values between groups in function ```TEgroup_stats.m``` can now be performed on raw or surrogate-corrected TE values (through parameter ```cfg.rawvalues```)
- Statistics script for large ```N``` ([```example_script_large_datasets.m```] [ex_script]) 

## 3.4 - 2015-10-07
This version of TRENTOOL was tested with fieldtrip-lite-20150928 and MATLAB version R2012B.

As of this version, we will no longer maintain parallel code in TRENTOOL to work with newer/future MATLAB versions (we are, however, happy to accept any code contributions).

### Changed
- ```InteractionDelayReconstruction_calculate.m``` performs a permutation test for the optimized delay ```u``` only; using this function is now much faster:
    - the function ```InteractionDelayReconstruction_analyze.m``` was replaced by ```TEfinddelay``` and ```TEfindmaxte```, called internally by IDR_calculate.m
    - there is only one file saved to disk containing all necessary information (TE raw data and statistical results), we also cleaned up the output of IDR_calculate, fixes issue [#9]
    - IDR_plotting can be used on the new output of IDR_calculate.m (this change is not backwards compatible), fixes issue [#14]
    
### Fixed
- The number of permutations is now calculated correctly in ```TEchecknumperm.m``` (issue [#5])

### Added
- Comparison of conditions within a single subject was added to the group analysis work flow (```TEgroup_conditionstatssingle.m```)
- Console outputs are now handled by a dedicated function ```TEconsoleoutput.m```, that allows the user to set the verbosity of TRENTOOL outputs (this has no effect on functionality)
- ```TEgetACT.m``` lets you calculate the autocorrelation time before starting an actual TE analysis
- ```TEchecknumperm.m```

### Removed
- ```InteractionDelayReconstruction_analyze.m```
- ```TEgroup_calculate.m``` (replaced by ```TEgroup_prepare.m``` and ```TEgroup_stats.m```)
- ```TEconditionstatssingle.m``` (replaced by ```TEgroup_conditionstatssingle.m```)
- ```TECvalues.m``` (contains outdated estimators)


[ex_script]: https://github.com/trentool/TRENTOOL3/blob/master/example_script_large_datasets.m
[#16]: https://github.com/trentool/TRENTOOL3/issues/16


[#5]:  https://github.com/trentool/TRENTOOL3/issues/5
[#9]:  https://github.com/trentool/TRENTOOL3/issues/9
[#14]: https://github.com/trentool/TRENTOOL3/issues/14


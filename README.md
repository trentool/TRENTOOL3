# TRENTOOL3
The **Tr**ansfer **En**tropy **Too**lbox, an Open-Source MATLAB toolbox for transfer entropy estimation, Reference: [Lindner et al. (2011)](http://www.biomedcentral.com/1471-2202/12/119).

TRENTOOL is an open-source MATLAB toolbox that allows the user to easily handle the considerable complexity of transfer entropy (TE) estimation from time series. For the use with neural data TRENTOOL seamlessly
integrates with the popular [FieldTrip toolbox](http://www.fieldtriptoolbox.org/).

TRENTOOL provides the following features:
* Transfer entropy estimation 
* Reconstruction of interaction delays 
* Time-resolved transfer entropy estimation
* Group analysis and statistics
* Partial correction for multivariate interactions

## Implementation
TRENTOOL is implemented in Mathworks MATLAB (The MathWorks Inc., Natick, MA, 2008) with some functions written in NVIDIA CUDA C/C++ code (NVIDIA Corporation, 2013). TRENTOOL also makes use of the MATLAB toolboxes [FieldTrip](http://www.fieldtriptoolbox.org/) and [TSTOOL](http://www.dpi.physik.uni-goettingen.de/tstool/). The user interacts with TRENTOOL through MATLAB scripts (.m-files).

## Releases
For a quick start download the latest release v3.4.2 [here](https://github.com/trentool/TRENTOOL3/releases).

All changes are documented in the [change log](https://github.com/trentool/TRENTOOL3/blob/master/CHANGELOG.md).

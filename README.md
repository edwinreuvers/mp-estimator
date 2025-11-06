# Muscle Property Estimator

This repository contains a toolbox (+ example scripts) to estimate muscle properties from quick-release, step-ramp and isometric experiments.

![Python Version](https://img.shields.io/badge/python-3.x-blue) ![License](https://img.shields.io/badge/license-MIT-green)

------------------------------------------------------------------------

## Table of Contents

-   [Installation](#installation)
-   [Dependencies](#dependencies)
-   [Usage](#usage)
-   [Citation](#citation)
-   [Contact / Support](#contact--support)

## Installation

1.  **Download the folder** `mp_estimator`.

2.  **Use directly** You can use the functions from within the same folder by importing them in Python:

``` python
from mp_estimator import estimate, loaddata
```

3.  **Install in your virtual environment** To make the package available anywhere in your virtual environment, navigate to the folder containing `setup.py` and run:

``` bash
pip install -e .
```

-   The `-e` flag means “editable,” so any changes you make to the code will be reflected immediately.
-   After this, you can import your functions anywhere in that environment:

``` python
from mp_estimator import estimate, loaddata
```

## Dependencies

This toolbox requires the following Python packages:

- [NumPy](https://numpy.org/)  
- [SciPy](https://www.scipy.org/)  
- [pandas](https://pandas.pydata.org/)

## Usage

The open-source toolbox automates muscle parameter estimation using data from quick-release, step-ramp, and isometric experiments. Our associated manuscript demonstrates that accuracy improves when correcting for CE shortening during quick-release. However, this requires both quick-release and step-ramp data. Recognising that some users may only have one type of data, we designed the toolbox with flexibility. The toolbox offers the following options:

1.  **Estimate CE, SEE and PEE force-length parameters only**
    -   Requires *only* quick-release data.
    -   Does **not** correct for CE shortening due to the quick-release.
2.  **Estimate CE force-velocity parameters only**
    -   Requires *only* step-ramp data.
    -   Since force-length parameter values are unknown, the following assumptions are made:
        a)  CE velocity during ramp equals MTC velocity
        b)  CE force equals SEE force
        c)  CE length during ramp equals CE optimum length
3.  **Estimate force-length and CE force-velocity parameters using the 'traditional method'**.
    -   Requires *both* quick-release and step-ramp data.
    -   Does **not** correct for CE shortening due to the quick-release.
4.  **Estimate force-length and CE force-velocity parameters using the 'improved method'** (recommended)
    -   Requires *both* quick-release and step-ramp data.
    -   **Does** correct for CE shortening due to the quick-release.
5.  **Estimate (de-)activation dynamics parameters**
    -   Any data type may be used, but isometric contraction data is recommended to minimise the influence of contraction dynamics (see manuscript).
    -   Can only be performed after estimating force-length and force-velocity parameters (via option 3 or 4).

The following sections explain how to use each option. Sample scripts are available for all options to further assist with implementation.

### CE, SEE and PEE force-length parameter value estimation

First, extract the relevant data from the quick-release experiments by setting the configurations below:

``` python
optsQR = {
    'dataDir':          'C:/xxx/xxx/xxx/Data/QR/',
    'iCols':            (0,1,3),
    'idxQRmin':         [(0,20), (1,11), (60,85)],
    'idxQRpre':         'auto',
    'idxQRpst':         'auto',
    'nQRsamp':          'auto',
    'dispFig':          True,
    }
```

-   **dataDir:** Path to your quick-release data files. If omitted, a file selector dialog will appear.
-   **iCols:** Tuple indicating the columns representing 1) time; 2 MTC length and 3) SEE force, respectively. If omitted, a column selector dialog will appear.
-   **idxQRmin:** List of start and stop indices for the period with minimum SEE force. Provide ranges for each quick-release experiment or use 'auto' for automatic detection.
-   **idxQRpre:** List of start and stop indices for SEE force before the quick-release. Provide ranges for each quick-release experiment or use 'auto' for automatic detection.
-   **idxQRpst:** List of start and stop indices for SEE force after the quick-release. Provide ranges for each quick-release experiment or use 'auto' for automatic detection.
-   **nQRsamp:** Number of samples to average before and after the quick-release. Default is 5 if 'auto' or unspecified.
-   **dispFig:** When True, plots will be displayed for manual verification and adjustment of selected data segments (i.e., idxQRmin, idxQRpre and idxQRpst).

Second, call the function 'loaddata.qr' to select the quick-release data:

``` python
dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr(optsQR)
```

-   **dataQR:** List of dictionaries containing the extracted quick-release data.
-   **idxQRmin, idxQRpre and idxQRpst:** Returned index ranges for your reference or future use.

Third, obtain the 'default' muscle parameter (such as width of the CE force-length relationship etc.):

``` python
defpar = loaddata.get_muspar()
```

-   **defpar:** Dictionary containing the 'default' muscle parameter values.

Fourth, estimate the CE, PEE and SEE and force-length parameter values by passing your data and initial parameters to 'estimate.fl':

``` python
estpar,dataQRout = estimate.fl(dataQR,defpar)
```

-   **estpar:** Dictionary containing estimated muscle parameter values.
-   **dataQRout:** List of dictionaries for each quick-release trial, with keys:
    -   *lseeQRpre:* Estimated SEE length before quick-release
    -   *lseeQRpst:* Estimated SEE length after quick-release
    -   *fseeQRpre:* SEE force before quick-release
    -   *fseeQRpst:* SEE force after quick-release
    -   *lpeeQR:* (Estimated) PEE length of the data
    -   *fpeeQR:* PEE force of the data
    -   *lmtcQRpre:* MTC length before quick-release

### CE force-velocity parameter value estimation

First, extract the relevant data from the step-ramp experiments by setting the configurations below:

``` python
optsSR = {
    'dataDir':          'C:/xxx/xxx/xxx/Data/SR/',
    'iCols':            (0,1,3),
    'idxSRcon':         'auto',
    'nSRsamp':          10,
    'dispFig':          False,
    }
```

-   **dataDir:** Path to your step-ramp data files. If omitted, a file selector dialog will appear.
-   **iCols:** Tuple indicating the columns representing 1) time; 2 MTC length and 3) SEE force, respectively. If omitted, a column selector dialog will appear.
-   **idxSRcon:** List of start and stop indices where SEE force is most ‘constant’. If set to 'auto', indices are automatically detected.
-   **nSRsamp:** Number of data points to average during the constant force phase. Default is 10 if 'auto' or unspecified.
-   **dispFig:** When True, plots will be displayed for manual verification and adjustment of selected data segments (i.e., idxSRcon).

Second, call the function 'loaddata.sr' to select the step-ramp data:

``` python
dataSR,idxSRcon = loaddata.sr(optsSR)
```

-   **dataSR:** List of dictionaries containing the extracted step-ramp data.
-   **idxSRcon:** Returned index ranges for your reference or future use.

Third, obtain the 'default' muscle parameter (such as width of the CE force-length relationship etc.):

``` python
defpar = loaddata.get_muspar()
```

-   **defpar:** Dictionary containing the 'default' muscle parameter values.

Fourth, estimate the CE force-velocity parameter values by passing your data and initial parameters to 'estimate.fv':

``` python
estpar,dataSRout = estimate.fv(dataSR,defpar)
```

-   **estpar:** Dictionary containing the estimated muscle parameters.
    -   If only CE force-velocity parameters are estimated (option 2), $a$, $b$, and $F_{CE}^{max}$ are estimated on based step-ramp data.
    -   If combined with estimation of force-length parameters (option 3 & 4), $a$, $b$ are estimated on based step-ramp data (while $F_{CE}^{max}$ is estimated based on quick-release data).
-   **dataQRout:** List of dictionaries for each quick-release trial, with keys:
    -   *vceSR:* Estimated CE shortening velocity.
    -   *fceSR:* Estimated CE force.
    -   *lcerelSR:* Relative CE length.

### Traditional method

Running the traditional method requires first estimating the CE, SEE and PEE force-length parameter values and then CE force-velocity parameter values. All inputs and outputs of the function are described above.

``` python
#%% Get default muscle parameters
defpar = loaddata.get_muspar()

#%% Settings and load data
optsQR = {
    # Include your options here
    }
optsSR = {
    #Include your options here
    }
dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr(optsQR)
dataSR,idxSRcon = loaddata.sr(optsSR)

#%% Run uncorrected method
estpar,dataQRout = estimate.fl(dataQR,defpar)
estpar,dataSRout = estimate.fv(dataSR,defpar,estpar)
```

### Improved method

To run the improved method, call a separate function ‘estimate.im’. This function iteratively call 'estimate.fl' and 'estimate.fv' until the change in all estimated parameter values are below 0.1%. All inputs and outputs of the function are described above.

``` python
#%% Get default muscle parameters
defpar = loaddata.get_muspar()

#%% Settings and load data
optsQR = {
    # Include your options here
    }
optsSR = {
    #Include your options here
    }
dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr(optsQR)
dataSR,idxSRcon = loaddata.sr(optsSR)

#%% Run improved method
estpar,dataQRout,dataSRout = estimate.im(dataQR,dataSR,defpar)
```

### Excitation dynamics parameter value estimation

To estimate the time constants of the calcium dynamics, one has to first estimate the parameter values of the CE, PEE and SEE force-length relationships as well as the parameters of the CE force-velocity relationship. After this the time constant of the calcium dynamics can be estimated. First, extract the relevant data (from the quick-release experiments) by setting the configurations below:

``` python
optsISOM = {
    'dataDir':          dataDir+'ISOM',
    'iCols':            (0,1,3,2),
    'durStimOffset':    0.1,
    'dispFig':          False,
    }
```

-   **dataDir:** Path to your (isometric) data files. If omitted, a file selector dialog will appear.
-   **iCols:** Tuple indicating the columns representing 1) time; 2 MTC length and 3) SEE force, respectively. If omitted, a column selector dialog will appear.
-   **durStimOffset:** Time (in seconds) to include after stimulation offset. If not specified, the default is 0.1 s.
-   **dispFig:** When True, plots will be displayed for manual verification of selected data segments (i.e., idxSRcon).

Second, call the function 'SelISOMdata' to select the (isometric) data:

``` python
dataISOM, idxSEL = SelISOMdata(optsISOM)
```

-   **dataISOM:** List of dictionaries containing the extracted (isometric) data.
-   **idxSEL:** List of start and stop indices of selected interval of the data.

Third, estimate the excitation dynamics parameter values by passing your data and initial parameters to 'getACTparms':

``` python
estpar,dataACTout = estimate.act(dataISOM,defpar,estpar)  
```

-   **estpar:** Dictionary containing the estimated muscle parameters.
-   **dataACTout:** List of dictionaries for each (isometric) experment, with keys::
    -   *time:* Time-axis.
    -   *lmtc:* MTC length.
    -   *fseeData:* Experimental SEE force.
    -   *fseeMdl:* Model-predicted SEE force.
    -   *tStim:* List containing stimulation onset and offset times.

## Citation

If you use this toolbox in your research or projects, please cite the [paper](https://doi.org/10.1101/2025.09.29.678508):

``` bibtex
@article{reuvers_2025_mpe,
  title={Accuracy of experimentally estimated muscle properties: Evaluation and improvement using a newly developed toolbox},
  author={Reuvers, Edwin D. H. M. and Kistemaker, Dinant A.},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.09.29.678508}
}
```

## Contact / Support

If you have questions, run into issues, or need help using the toolbox, you can reach out to me:

-   **GitHub Issues:** <https://github.com/edwinreuvers/mp-estimator/issues>
-   **Bluesky:** <https://bsky.app/profile/edwinreuvers.bsky.social>

# pafpcluster

**MATLAB routines for estimating number of molecules from a localization data of 
super-resolution localization microscopy**

## Prerequisites
1. The fluorescence kinetic rates (photobleaching, dark state entrance, dark state
exit) of the PA-FP used for experiment should be estimated, for example, through in vitro measurement.   

2. All the PA-FP molecules should be photoactivated more or less at a constant rate within
a given experimental time. 

## What it does
* **optimalsptcluster**: Iterative solver for both the optimal dark time and the number
of molecules  

* **fastsptcluster**: Fast Hoshen-Kopelman clustering

* **optmataucsolve2**: Approximate solver for optimal fluorescence dark time to correct
blinking-induced molecular overcounting. The kinetic rates, experimental time, and 
the number of molecules should be provided as input.  

* **optmtaucsolve2_dendra**: Approximate optimal dark time solver using the rates measured in vitro with the microscope in Bustamante lab. 

* **sampledata.mat**: Sample data

* **samplecode**: Sample script

## Usage
A sample localization data and a script is provided. Simply run 'samplecode' in
matlab environment. 

## References
1. J. Y. Shin, J. Lopez-Garrido, S-. H. Lee, C. Diaz-Celis, T. Fleming, 
C. Bustamante and K. Pogliano
"Visualization and functional dissection of coaxial paired SpoIIIE channels across the sporulation septum,"
_eLIFE_ **e06474** (2015).

2. S-.H. Lee, J. Y. Shin, A. Lee and C. Bustamante, 
"Counting single photoactivatable fluorescent molecules by 
photoactivated localization microscopy (PALM),"
_Proc. Natl. Acad. Sci._ **109**, 17436-17441 (2012).
 
3. A. Al-Futaisi and T. W.. Patzek, 
"Extension of Hoshen-Kopelman algorithm to non-lattice environments," 
_Physica A_ **321**, 665-678 (2003).




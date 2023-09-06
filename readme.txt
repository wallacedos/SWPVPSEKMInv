The code package is expected to solve the mode-losses and mode-aliases of
the inversion for the multimodal surface wave dispersion curves.
 
Idea and operation:
The pattern search (PS) is used to invert the reliable segment of the 
fundamental-mode surface wave phase velocities for the first stage. For 
the second stage, the inverted result of the first stage is set as the 
initial model, the PS with embedded Kuhn-Munkres (PSEKM) algorithm is 
adopted for inverting the observed phase velocities of all modes. And for 
each frequency, a weighted bipartite graph is established between the 
observed values with no-explicitly-specified-mode-order (NESMO) and 
predicted values of the model m during the inversion, then the maximum 
match is determined by the Kuhn-Munkres algorithm for calculating the 
minimum distance between the observed and predicted data sets. The 
mode-order information of the observed phase velocities with NESMO would 
be dynamically evaluated for each model m occurred in the inversion process.

Note:
Before using the inversion workflow, we strongly recommend that you read 
our article and the tutorial. The article and tutorial are placed in 
"main-folder"\Doc. If you publish an article by using our inversion programs, 
please cite our article. If there are any questions about the code, please 
contact with the first author Yan Yingwei by email "wallace2012y@outlook.com".

The code has been debugged on the the matlab R2019a. A higher or lower 
version of matlab may also make the code run successfully.

"main-folder"\Example consists of three sample scripts.

roadBed1.m: reproduction of the inversion of roadbed 1 of our article 1.

DispersionMeasurement_GVDA.m: Measuring the multimodal dispersion of the
distributed-acoustic-sensing (DAS) data collected at Long Line I at 
Garner Vallry Downhole Array (GVDA), California. We utlize the cylindrical-wave
phase to measure the multimodal dispersion from the CCFs of the DAS data.

DASDCDataInversion_GVDA.m: Inversion of the multimodal dispersion curves of the 
DAS data. We adopt the modern inversion workflow in our article. 

DispersionMeasurement_GVDA.m and DASDCDataInversion_GVDA.m reproduce the result
of our article 2.

Our article:
Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
the multimodal surface wave dispersion curves: Staging strategy and Pattern 
search with embedded Kuhn-Munkres algorithm, Geophysical Journal
International,231(01), 47-71,  
https://doi.org/10.1093/gji/ggac178. 

Yan， Y., Chen, X., Li, J., Guan, J., Xi, C., Liu, H. 2023. Inversion of
multimodal dispersion curves from distributed acoustic sensing measurements
for subsurface imaging: A field case of Garner Valley, California, Journal of
Applied Geophysics, 214, 105070,
https://doi.org/10.1016/j.jappgeo.2023.105070 

Author(s): Yan Yingwei
Email:     wallace2012y@outlook.com
Copyright: 2022-2025 
Revision:  1.0  Date: 9/6/2023

Department of Earth and Space Sciences, Southern University of Science 
and Technology (SUSTech).
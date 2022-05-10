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
our article. The article is placed in "SWPVPSEKMInv\Doc". If you publish
an article by using our inversion programs, please cite our article. If
there are any questions about the code, please contact with the first 
author Yan Yingwei by email "wallace2012y@outlook.com".
The code has been debugged on the the matlab R2019a. A higher or lower 
version of matlab may also make the code run successfully.

Our article:
Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of 
the multimodal surface wave dispersion curves: Staging strategy and Pattern 
search with embedded Kuhn-Munkres algorithm, Geophysical Journal
International,ggac178, 
https://doi.org/10.1093/gji/ggac178. 

Author(s): Yan Yingwei
Email:     wallace2012y@outlook.com
Copyright: 2022-2025 
Revision:  1.0  Date: 5/10/2022

Department of Earth and Space Sciences, Southern University of Science 
and Technology (SUSTech).
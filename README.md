# WavePart ![image](https://user-images.githubusercontent.com/48567321/126871118-68ae6e82-15fe-48e5-ac1e-2d1d0b673a9b.png)


Set of Matlab(r) functions for the partition of directional wave spectra to its wind and different swell components. The partitions are initially identified using a watershed defining algorithm  and then are modified following mostly the method described in Hanson and Phillips (2001).

Authors:  
  Douglas Cahl and George Voulgaris  
  School of the Earth, Ocean and Environment  
  University of South Carolina  
  Columbia, SC, 29205, USA.  
  
Code Updates:
  -  1/22/2020 - waveparams.m - the mean freq. (fm) estimate was incorrect; it has been corrected.
  -  1/22/2020 - partition.m was updated to catch cases when a flat spectrum is given as input.
  -  1/22/2020 - master_partition.m was updated so that it calls waveparams.m with the correct number of outputs.
  -  1/21/2020 - waveparams.m - the Hrms estimate was incorrect; also the mean direction now is given in -180 to +180 deg range.
  
Citation:  
   -  Douglas, C. and Voulgaris, G., 2019, WavePART: MATLAB(r) software for the partition of directional ocean wave spectra.: Zenodo, doi:10.5281/zenodo.2638500. 

Relevant References:  
   -  J.L. Hanson and O.M. Philips, 2001. Automated Analysis of Ocean Surface Directional  Wave Spectra. Journal of Oceanic and Atmospheric Technology, 18, 278-293.   
   -  J. Portilla, F.J. Ocampo-Torres, and J. Monbaliu, 2009. Spectral Partitioning and Identification of Wind Sea and Swell.  Journal of Oceanic and Atmospheric Technology, 26, 107-121. DOI: 10.1175/2008JTECHO609.1   
   -  E. Cheynet, 2019. Pcolor in Polar Coordinates (https://www.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates), MATLAB Central File Exchange. Retrieved March 16, 2019.  

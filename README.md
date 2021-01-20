# MORB-CO2-vesicles

This repository contains code for reproducing stereological analysis of images and x-ray microtomography scans in  <em>Jones, M., Soule, S., Liao, Y., Brodsky, H., Le Roux, V., Klein, F., Quantitative vesicle analyses and total CO<sub>2</sub> reconstruction in mid-ocean ridge basalts</em> (in revision at the Journal of Volcanology and Geothermal Research).

## Authors

Meghan Jones<sup>a</sup>, Adam Soule<sup>b</sup>, Yang Liao<sup>b</sup>, Harry Brodsky<sup>c</sup>, Veronique Le Roux<sup>b</sup>, Frieder Klein<sup>b</sup>

<sup>a</sup>Massachusetts Institute of Technology/Woods Hole Oceanographic Institution Joint Program in Oceanography, 360 Woods Hole Road, 02543, Woods Hole, MA, USA

<sup>b</sup>Department of Geology and Geophysics, Woods Hole Oceanographic Institution, 02543, Woods Hole, MA, USA

<sup>c</sup>Department of Mechanical and Industrial Engineering, Northeastern University, 02115, Boston, MA, USA

## Key Points

* Stereological methods closely reproduce 3D MORB vesicle size distributions and vesicularities measured by x-ray micro-tomography
* Revised methods for calculating total CO<sub>2</sub> concentrations in MORB provide more accurate estimates
* New total CO<sub>2</sub> estimates have implications for estimating mantle carbon concentrations and ridge CO<sub>2</sub> flux 

## Abstract

Vesicle textures in submarine lavas have been used to calculate total (pre-eruption) volatile concentrations in mid-ocean ridge basalts (MORB), which provide constraints on upper mantle volatile concentrations and global mid-ocean ridge CO<sub>2</sub> flux. In this study, we evaluate vesicle size distributions (VSDs) and volatile concentrations in a suite of 20 MORB samples that span the range of vesicularities and vesicle number densities observed in MORB globally. We provide recommended best practices for quantifying vesicularity, vesicle number densities, and VSDs based on synthetic vesicle populations and comparisons between traditional 2D methods and x-ray computed micro-tomography results. For 2D measurements, we recommend analyzing multiple polished fragments with a cumulative area >100 times the area of the largest observed vesicle and including >200 vesicles in stereological VSD reconstructions. For 3D measurements, we recommend analyzing sample volumes >0.01 cm<sup>3</sup> at resolutions <2.0 μm/pixel for low vesicularity MORB (i.e., <4 vol.%) and sample volumes >0.1 cm<sup>3</sup> with resolutions <5 μm/pixel for higher vesicularity samples. Our validation of vesicularity measurements allows reconstructions of total CO<sub>2</sub> concentrations in MORB using dissolved volatile concentrations, vesicularities, and equations of state. We assess approaches for estimating the exsolved CO<sub>2</sub> concentration in MORB vesicles and find that CO<sub>2(g)</sub> density is ~40% lower than previously suggested, likely due to melt contraction during quenching. Based on these results, we recommend using sample eruption pressures, magmatic temperatures, and an equation of state that accounts for non-ideality at high temperatures to calculate exsolved CO<sub>2</sub> when independent constraints from Raman spectroscopy or laser ablation are unavailable. Our results suggest that some previous studies may have overestimated MORB volatile concentrations by as much as 50%, with the greatest differences in samples with the highest vesicularities. These new results imply lower CO<sub>2</sub>/Ba of undegassed, enriched-MORB and lower integrated global ridge CO<sub>2</sub> flux than previously inferred.

## Status

This paper is in revision at the Journal of Volcanology and Geothermal Research. Additional code to reproduce the synthetic vesicle and Raman spectroscopy analyses and comparisons between microCT resolutions will be added to this repository prior to resubmission. Comments, questions, or suggestions are appreciated through github issues or via e-mail to Meghan Jones (meghanj [at] alum.mit.edu)

## Acknowledgments and support 

We are thankful to the captain, crew, vehicle teams, and science participants of the R/V Thompson VISIONS’11 cruise, R/V Western Flyer Northern Expeditions cruises, and R/V Atlantis AT33-03 cruise for assistance in collecting the samples used in this study. We thank T. Grove, D. Lizarralde, M. Kurz, T. Perron, M. Manga, and D. Wanless for insightful comments on earlier versions of the manuscript. We thank R. Bodnar, G. Gaetani, H. Lamadrid, A. Pamukcu, and T. Shea for helpful conversations. M. Jones was supported by the Department of Defense (DoD) through the National Defense Science & Engineering Graduate Fellowship (NDSEG) Program. This work was supported by an ExxonMobil grant and NSF grants OCE-1333492, OCE-1259218, and OCE-1260578.

## Contents

* Images and derived data (from ImageJ and Bruker Micro-CT software) for a subset of the images/scans used in the manuscript.
* Code to load data from images and x-ray micro-tomography scans and process the data using the stereological methods from Cheng and Lemlich (1983), Sahagian and Proussevitch (1998), and Saltikov (1967). 
* Jupyter notebook showing an example of stereological data correction for sample AX13-RC04.

## Citations for equations used in code

Cheng, H.C., Lemlich, R., 1983. Errors in the measurement of bubble size distribution in foam. Ind. Eng. Chem. Fundam. 22, 105–109. https://doi.org/10.1021/i100009a018
    
Sahagian, D.L., Proussevitch, A.A., 1998. 3D particle size distributions from 2D observations: stereology for natural applications. J. Volcanol. Geotherm. Res. 84, 173–196. https://doi.org/10.1016/S0377-0273(98)00043-2

Saltikov, S., 1967. The determination of the size distribution of particles in an opaque material from a measurement of the size distribution of their sections, in: Stereology: Proceeding of the Second International Congress for Stereology. pp. 163–173.

Pre-processing for the data included image processing in MATLAB, Bruker MicroCT Software, and FIJI:

Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch, T., Preibisch, S., Rueden, C., Saalfeld, S., Schmid, B., Tinevez, J.-Y., White, D.J., Hartenstein, V., Eliceiri, K., Tomancak, P., Cardona, A., 2012. Fiji: an open-source platform for biological-image analysis. Nat. Methods 9, 676–682. https://doi.org/10.1038/nmeth.2019
        
Early development of this code also referenced the FOAMS MATLAB program:

Shea, T., Houghton, B.F., Gurioli, L., Cashman, K.V., Hammer, J.E., Hobden, B.J., 2010. Textural studies of vesicles in volcanic rocks: An integrated methodology. J. Volcanol. Geotherm. Res. 190, 271–289. https://doi.org/10.1016/j.jvolgeores.2009.12.003

Additional citations are referenced in the full manuscript. 

### Dependencies
  - python=3.8 
  - numpy
  - matplotlib
  - pandas
  - jupyter

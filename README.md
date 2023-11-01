
  # Check out the data challenge [website](https://2024b-dc-swotdenoising.readthedocs.io) for more infos !

<p align="center">
  <img src="figures/dc_2024b_SWOTdenoising_banner.jpg" alt="Alt Text" width="900"/>
</p>

# SWOT denoising data challenge
 
This repository contains codes and sample notebooks for downloading and processing the 2024b SWOT denoising data challenge.
Note that this data challenge is part of the extended effort to improve oceanographic algorithms using data challenges: [ocean-data-challenges](https://ocean-data-challenges.github.io/index.html).

So far, the github page visits amount to: 

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Focean-data-challenges%2F2024b_DC_SWOTdenoising&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising)


 
# 1. Context 


## The SWOT mission  

The Surface Water and Ocean Topography (SWOT) satellite mission is an innovative project designed to measure various aspects of Earth's water bodies, including the height of surface water, from oceans to land.

SWOT aims to measure the height of almost all surface water on Earth. It will provide data on oceans, land, and freshwater resources. This mission covers the entire water supply-demand chain, making it unique in its capability to monitor Earth's surface water entirely. SWOT offers very high-definition data, this level of detail enables scientists to monitor sea level changes near coasts, observe dynamic features in the open ocean, and track factors like ocean circulation and the transport of heat, energy, oxygen, and nutrients.
  
The mission involves extensive calibration and validation efforts. Field campaigns have been planned worldwide to compare SWOT observations with ground-based measurements. Scientists from various regions are coming together to contribute to the calibration phase. Although the original requirement was to deliver validated datasets a year after SWOT's launch, NASA and CNES are exploring a new approach to distribute data earlier to the science community. This approach encourages active participation from the scientific community, through the SWOT Science Team, in the validation phase.
 
 

## The intrumental KaRIn noise
 

The instrumental noise in KaRIn (Ku-band Radar Interferometer), the key instrument onboard the SWOT satellite, is a critical factor to consider and mitigate in the processing of SWOT data. KaRIn is designed to measure the height of Earth's water surfaces with high precision. To achieve the mission's goals, it needs to be extremely sensitive to variations in sea surface height (SSH). However, like any instrument, KaRIn is subject to various sources of noise that can affect the accuracy of its measurements.
This instrumental noise can introduce errors or uncertainty into the measurements, which need to be corrected or removed to ensure the data's accuracy.

In regions with low noise, the impact on the data may be minimal. However, in areas with more significant noise, the presence of instrumental noise can lead to inaccuracies in the derived measurements of SSH.For instance, surface wave height can be a factor. High surface wave conditions can make it challenging to obtain precise measurements of water height. In such conditions, the instrumental noise may become more pronounced, making it difficult to remove it effectively. Also, the quality of first and second derivatives of SSH data along SWOT tracks is crucial for understanding ocean dynamics, such as currents and eddies. Errors introduced by instrumental noise can impact the derivatives and hinder the accuracy of these derived oceanographic parameters.


The Unet method used for noise removal is a valuable tool for mitigating instrumental noise in SWOT data. However, as with any data processing technique, it may have limitations, particularly in regions with extreme conditions (that were not seen in the Unet's training) such as high surface wave heights. When the Unet correction fails to reduce noise in some areas, it indicates the need for more robust methods or improvements in the processing techniques.


To address this challenge, researchers and scientists working with SWOT data need to continuously refine and improve the noise removal and data processing methods. This includes developing more sophisticated algorithms, considering the unique conditions in different regions, and exploring ways to optimize the Unet method or develop new approaches for noise removal. 

This is precisely why initiatives like the present data challenge are crucial in the context of the SWOT mission. The challenges related to instrumental noise in KaRIn data, as well as the difficulties encountered in processing SWOT data in regions with unique conditions highlight the need for collaborative efforts and innovative solutions.

By organizing this data challenge, experts and researchers can come together to share their expertise, test various noise removal and data processing methods, and collectively work towards more accurate and robust results. This collaborative approach fosters innovation and encourages the development of improved techniques for handling instrumental noise and deriving valuable information from SWOT data.

# 2. The data challenge setup 
 

## The SWOT denoising regions


<embed>  
<center>
<div class="image-container"> 
<map name="map_example1"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_acc.html" target="_blank" alt="ACC" shape=poly coords="335,177, 335,165, 348,165, 348,177"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_gulfstream.html" target="_blank" alt="Gulf Stream" shape=poly coords="157,95, 157,83, 170,83, 170,95"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_med.html" target="_blank" alt="Mediterranean Sea" shape=poly coords="203,90, 203,80, 215,80, 215,90"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_cali.html" target="_blank" alt="Californian Current" shape=poly coords="92,95, 92,83, 105,83, 105,95">  <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_mada.html" target="_blank" alt="Madagascar shelf" shape=poly coords="243,132, 243,121, 255,121, 255,132">  <img src="https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising/assets/33433820/b62394ab-cc2b-488e-82a7-421474e48b84" title="1d orbit" alt="image map example" width=400 height=250 usemap="#map_example1"></map>
<map name="map_example2"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_acc.html" target="_blank" alt="ACC" shape=poly coords="335,177, 335,165, 348,165, 348,177"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_gulfstream.html" target="_blank" alt="Gulf Stream" shape=poly coords="157,95, 157,83, 170,83, 170,95"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_med.html" target="_blank" alt="Mediterranean Sea" shape=poly coords="203,90, 203,80, 215,80, 215,90"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_cali.html" target="_blank" alt="Californian Current" shape=poly coords="92,95, 92,83, 105,83, 105,95"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_mada.html" target="_blank" alt="Madagascar shelf" shape=poly coords="243,132, 243,121, 255,121, 255,132"> <img src="https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising/assets/33433820/7cef9610-3380-4c0c-abde-b706b1b6518f" title="21d orbit" alt="image map example" width=400 height=250 usemap="#map_example2"></map> 
</div> 
</center>
</embed>


* Gulf Stream: [1 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_gulfstream.html) and [21 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_gulfstream.html)
    - Area: [55°W, 45°W, 30°N, 40°N]
    - Dynamical specificities: High variability region with mixed geostrophic and ageostrophic dynamics.
    
* Mediterranean Sea: [1 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_med.html) and [21 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_med.html)
    - Area: [2°W, 8°E, 36°N, 44°N]
    - Dynamical specificities: A quasi-closed basin with strong ageostrophic dynamics and vertical shear.

* Atlantic circumpolar [1 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_acc.html) and [21 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_acc.html)
    - Area: [152°E, 162°E, 64°S, 54°S]
    - Dynamical specificities: Strong Surface Wave Height (SWH) region.
    
* Californian Current [1 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_cali.html) and [21 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_cali.html)
    - Area: [130°W, 120°W, 30°N, 40°N]
    - Dynamical specificities: Coastal high variability region with internal waves.
    
* Madagascar shelf [1 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_1dswotdenoising/1dswotdenoising_mada.html) and [21 day orbit](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/3_21dswotdenoising/21dswotdenoising_mada.html)
    - Area: [44°W, 54°W, 12°S, 2°S]
    - Dynamical specificities: Shelf with shallow topography and inaccurate MSS representation.


## Input data

**SWOT file structure**

<p align="center">
  <img src="figures/illust_SWOT_netcdffile.png" alt="Alt Text" width="500"/>
</p>

 
## Evaluation

- [Visualisation](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/4_metrics_det/metrics_1-visualisation.html)

- [Physical pdf](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/4_metrics_det/metrics_2-physicalpdf.html)

- [Spectrum](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/4_metrics_det/metrics_3-spectrum.html)

- [Energy loss](https://2024b-dc-swotdenoising.readthedocs.io/en/latest/4_metrics_det/metrics_4-energyloss.html)


## Denoising techniques

 
 

  

# 3. Get started
 

## Installation
:computer: _**How to get started ?**_

Clone the data challenge repo: 
```
git clone https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising.git
```
or using SSH: 
```
git clone git@github.com:ocean-data-challenges/2024b_DC_SWOTdenoising.git
```

create the data challenge conda environment, named env-dc-swot-filtering, by running the following command:
```
conda env create --file=dc_environment.yml 
```
and activate it with:

```
conda activate env-dc-swot-denoising
```
then add it to the available kernels for jupyter to see: 
```
ipython kernel install --name "env-dc-swot-denoising" --user
```
finally, select the "env-dc-woc-esa" kernel in your notebook with Kernel > Change Kernel.

You're now good to go ! 


## Download the data

The data are hosted and can be accessed on the MEOM server opendap [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/catalog.html). You can store the downloaded dataset in the swot_data/ directory. 

To download the SWOT dataset locally, you will need approximately 4 GB of disk space.  

 

  
  

# Acknowledgement


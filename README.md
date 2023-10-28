
  # Check out the data challenge [website](https://2024b-dc-swotdenoising.readthedocs.io) for more infos !

<p align="center">
  <img src="figures/dc_2024b_SWOTdenoising_banner.jpg" alt="Alt Text" width="900"/>
</p>

# SWOT denoising data challenge
 
This repository contains codes and sample notebooks for downloading and processing the 2024b SWOT denoising data challenge.
Note that this data challenge is part of the extended effort to improve oceanographic algorithms using data challenges: [ocean-data-challenges](https://ocean-data-challenges.github.io/index.html).

So far, the github page visits amount to: 

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Focean-data-challenges%2F2024_DC_WOC-ESA&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=PAGE+VIEWS&edge_flat=false)](https://github.com/ocean-data-challenges/2024_DC_WOC-ESA)



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

## Input data
 
## Evaluation

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

The data are hosted and can be accessed on the MEOM server opendap [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/catalog.html). The disk space needed to locally download the full dataset (for the reconstruction experiment, the independant evaluation and the comparison) is approximately 33Go. The comparison data is by far the heaviest with approximately 26Go. 

 

  
  

# Acknowledgement


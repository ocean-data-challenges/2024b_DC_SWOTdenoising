.. 2024_DC_WOC-ESA documentation master file, created by
   sphinx-quickstart on Fri Jul 21 14:53:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive. 
    
    
=============================
SWOT denoising data challenge
=============================

.. role:: raw-html(raw)
    :format: html

:raw-html:`<br />`

.. image:: ../../figures/dc_2024b_SWOTdenoising_banner.jpg
    :width: 600 
    :align: center

:raw-html:`<br />`

:raw-html:`<br />`

    
1. Context
==========
     

The SWOT mission 
----------------

The Surface Water and Ocean Topography (SWOT) satellite mission is an innovative project designed to measure various aspects of Earth's water bodies, including the height of surface water, from oceans to land.

SWOT aims to measure the height of almost all surface water on Earth. It will provide data on oceans, land, and freshwater resources. This mission covers the entire water supply-demand chain, making it unique in its capability to monitor Earth's surface water entirely. SWOT offers very high-definition data, this level of detail enables scientists to monitor sea level changes near coasts, observe dynamic features in the open ocean, and track factors like ocean circulation and the transport of heat, energy, oxygen, and nutrients.
  
The mission involves extensive calibration and validation efforts. Field campaigns have been planned worldwide to compare SWOT observations with ground-based measurements. Scientists from various regions are coming together to contribute to the calibration phase. Although the original requirement was to deliver validated datasets a year after SWOT's launch, NASA and CNES are exploring a new approach to distribute data earlier to the science community. This approach encourages active participation from the scientific community, through the SWOT Science Team, in the validation phase.
 
  

The intrumental KaRIn noise
---------------------------

The instrumental noise in KaRIn (Ku-band Radar Interferometer), the key instrument onboard the SWOT satellite, is a critical factor to consider and mitigate in the processing of SWOT data. KaRIn is designed to measure the height of Earth's water surfaces with high precision. To achieve the mission's goals, it needs to be extremely sensitive to variations in sea surface height (SSH). However, like any instrument, KaRIn is subject to various sources of noise that can affect the accuracy of its measurements.
This instrumental noise can introduce errors or uncertainty into the measurements, which need to be corrected or removed to ensure the data's accuracy.

In regions with low noise, the impact on the data may be minimal. However, in areas with more significant noise, the presence of instrumental noise can lead to inaccuracies in the derived measurements of SSH.For instance, surface wave height can be a factor. High surface wave conditions can make it challenging to obtain precise measurements of water height. In such conditions, the instrumental noise may become more pronounced, making it difficult to remove it effectively. Also, the quality of first and second derivatives of SSH data along SWOT tracks is crucial for understanding ocean dynamics, such as currents and eddies. Errors introduced by instrumental noise can impact the derivatives and hinder the accuracy of these derived oceanographic parameters.


The Unet method used for noise removal is a valuable tool for mitigating instrumental noise in SWOT data. However, as with any data processing technique, it may have limitations, particularly in regions with extreme conditions (that were not seen in the Unet's training) such as high surface wave heights. When the Unet correction fails to reduce noise in some areas, it indicates the need for more robust methods or improvements in the processing techniques.


To address this challenge, researchers and scientists working with SWOT data need to continuously refine and improve the noise removal and data processing methods. This includes developing more sophisticated algorithms, considering the unique conditions in different regions, and exploring ways to optimize the Unet method or develop new approaches for noise removal. 

This is precisely why initiatives like the present data challenge are crucial in the context of the SWOT mission. The challenges related to instrumental noise in KaRIn data, as well as the difficulties encountered in processing SWOT data in regions with unique conditions highlight the need for collaborative efforts and innovative solutions.

By organizing this data challenge, experts and researchers can come together to share their expertise, test various noise removal and data processing methods, and collectively work towards more accurate and robust results. This collaborative approach fosters innovation and encourages the development of improved techniques for handling instrumental noise and deriving valuable information from SWOT data.
 
 

2. The data challenge setup
===========================

The SWOT denoising regions
-------------------------- 

* `Mediterranean Sea <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_med.html>`_
    * Area: [5°W, 25°E, 35°N, 47°N]
    * Dynamical specificities: A quasi-closed basin with strong ageostrophic dynamics and vertical shear.  

* `Gulf Stream <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_gulfstream.html>`_ 
    * Area: [75°W, 45°W, 20°N, 50°N]
    * Dynamical specificities: High variability region with mixed geostrophic and ageostrophic dynamics.  
    
* `Atlantic circumpolar <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_circum.html>`_ 
    * Area: [14°E, 35°E, 45°S, 30°S]
    * Dynamical specificities: Strong Surface Wave Height (SWH) region.  
  
:raw-html:`<br />` 
    
.. raw:: html

    <embed> 
        <center>
        <div id="image_map"> <map name="map_example"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_circum.html" target="_blank" alt="ACC" shape=poly coords="310,310, 310,290, 340,290, 340,310"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_gulfstream.html" target="_blank" alt="Gulf Stream" shape=poly coords="160,150, 160,110, 210,110, 210,150"> <area href="https://2024b-dc-swotdenoising.readthedocs.io/en/latest/2_swotdenoising/swotdenoising_med.html" target="_blank" alt="Mediterranean Sea" shape=poly coords="300,155, 300,115, 345,115, 345,155">  <img src="https://github.com/ocean-data-challenges/2024_DC_WOC-ESA/assets/33433820/4be7d3ae-ea45-4b5c-b049-06ef985024c0" title="Gulf Stream" alt="image map example" width=600 height=350 usemap="#map_example"></map> </div> </center>

    </embed>

:raw-html:`<br />` 


Input data
----------

Evaluation
----------

Denoising techniques
--------------------

:raw-html:`<br />` 
    


3. Get started
============== 

Installation
------------
**How to Get Started?**

Clone the data challenge repo by using HTTPS:

.. code-block:: bash

    git clone https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising.git

Or clone it using SSH:

.. code-block:: bash

    git clone git@github.com:ocean-data-challenges/2024b_DC_SWOTdenoising.git

Create the data challenge conda environment, named ``env-dc-swot-filtering``, by running the following command:

.. code-block:: bash

    conda env create --file=dc_environment.yml

And activate it with:

.. code-block:: bash

    conda activate env-dc-swot-denoising

Then add it to the available kernels for Jupyter to see:

.. code-block:: bash

    ipython kernel install --name "env-dc-swot-denoising" --user

Finally, select the "env-dc-woc-esa" kernel in your notebook with Kernel > Change Kernel.

You're now good to go!



Download the data
----------------- 

The data is hosted and can be accessed on the MEOM server opendap `here <https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/catalog.html>`_. You can store the downloaded dataset in the swot_data/ directory. 

To download the SWOT dataset locally, you will need approximately 4 GB of disk space.  

 

:raw-html:`<br />` 

.. raw:: html

    <embed>    
        So far, the github page visits amount to: <br> <br> <a href="https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising"><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Focean-data-challenges%2F2024b_DC_SWOTdenoising&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false"/></a>
        
    </embed>
    
----------------- 

:raw-html:`<br />`
 
    
.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Get started 

   1_getstarted/getstarted_install.md
   1_getstarted/getstarted_data.md 
   1_getstarted/getstarted_eval.md

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: SWOT denoising

   2_swotdenoising/swotdenoising_details.md
   2_swotdenoising/swotdenoising_overalleval.md
   2_swotdenoising/swotdenoising_med.md
   2_swotdenoising/swotdenoising_gulfstream.md
   2_swotdenoising/swotdenoising_circum.md  
  
.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Metrics details

   3_metrics_det/metrics_1-visualisation.md
   3_metrics_det/metrics_2-energyloss.md
   3_metrics_det/metrics_3-physicalpdf.md 
    

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Scripts

   4_scripts/modules.rst

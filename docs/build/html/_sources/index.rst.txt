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
    :width: 1000 
    :align: center

:raw-html:`<br />`

:raw-html:`<br />`
 
     

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

The SWOT denoising data challenge
---------------------------------

To address this challenge, researchers and scientists working with SWOT data need to continuously refine and improve the noise removal and data processing methods. This includes developing more sophisticated algorithms, considering the unique conditions in different regions, and exploring ways to optimize the Unet method or develop new approaches for noise removal. 

This is precisely why initiatives like the present data challenge are crucial in the context of the SWOT mission. The challenges related to instrumental noise in KaRIn data, as well as the difficulties encountered in processing SWOT data in regions with unique conditions highlight the need for collaborative efforts and innovative solutions.

By organizing this data challenge, experts and researchers can come together to share their expertise, test various noise removal and data processing methods, and collectively work towards more accurate and robust results. This collaborative approach fosters innovation and encourages the development of improved techniques for handling instrumental noise and deriving valuable information from SWOT data.

To start with this data challenge you can check the `Participation <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/1_getstarted/index_getstarted.html>`_  section to learn how to: 

#. `Install <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/1_getstarted/getstarted_install.html>`_ the challenge, 
#. `Download the data <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/1_getstarted/getstarted_data.html>`_, 
#. `Evaluate your denoising techniques <https://2024b-dc-swotdenoising.readthedocs.io/en/latest/1_getstarted/getstarted_eval.html>`_. 
 
  

     
More data challenges   
--------------------

If you are interested in more data challenges relating to oceanographic data (global altimetric mapping, SWOT preprocessing techniques ...), you can visit the ocean-data-challenges website. 
  
  
    
.. raw:: html  


    <embed>  
        
        <br />
        
        <center><a  href="https://ocean-data-challenges.github.io"/> <img src="_static/odc_webpage.jpg" alt="Alt Text" width="500"/></a></center>
        
        <center><a  href="https://ocean-data-challenges.github.io" alt="Alt Text"/> ocean-data-challenges.github.io </a></center>
        
        <br /> 
        
        <br />
        

 

:raw-html:`<br />` 

.. raw:: html

    <embed>    
        So far, the github page visits amount to: <br> <br> <a href="https://github.com/ocean-data-challenges/2024b_DC_SWOTdenoising"><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Focean-data-challenges%2F2024b_DC_SWOTdenoising&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false"/></a>
        
    </embed>
    
----------------- 

:raw-html:`<br />`
 
    
.. toctree:: 
   :hidden:
   :maxdepth: 1
   :caption: Participate 

   1_getstarted/index_getstarted.md

.. toctree:: 
   :hidden:
   :maxdepth: 1
   :caption: Experimental Setup

   2_1dswotdenoising_setup/1dswotdenoising_details.md  

.. toctree:: 
   :hidden:
   :maxdepth: 1
   :caption: Regional evaluations
 
   2_1dswotdenoising_V03/index_eval_V03.md   

.. toctree:: 
   :hidden:
   :maxdepth: 1
   :caption: Regional evaluations
 
   2_1dswotdenoising_V102/index_eval_V102.md   
    
  
.. toctree::
   :hidden: 
   :maxdepth: 1
   :caption: Metrics details

   4_metrics_det/index_metrics.md 
    
.. toctree::  
    :caption: Contact us
    :hidden:
    :maxdepth: 0 
    
    contactus.md  

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Scripts

   5_scripts/modules.rst

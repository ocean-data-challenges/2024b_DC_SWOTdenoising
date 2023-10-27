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

    
     

SWOT data 
--------- 
 

:raw-html:`<br />`

Instrumental errors: KaRIn noise
--------------------------------
 
The data challenge setup
------------------------

The SWOT denoising regions:  

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

:raw-html:`<br />` 

.. raw:: html

    <embed>    
        So far, the github page visits amount to: <br> <br> <a href="https://github.com/ocean-data-challenges/2024_DC_WOC-ESA"><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Focean-data-challenges%2F2024_DC_WOC-ESA&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=PAGE+VIEWS&edge_flat=false"/></a> 
        
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

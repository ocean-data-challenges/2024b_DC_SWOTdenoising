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

* `Mediterranean Sea <https://2024-dc-woc-esa.readthedocs.io/en/latest/2_dc_medsea/dc_medsea_details.html>`_
    * Area: [5°W, 25°E, 35°N, 47°N]
    * Dynamical specificities: A quasi-closed basin with strong ageostrophic dynamics and vertical shear. 
    * WOC products: Merged SSH/SST currents, Drifter data-driven currents, BFN-QG geostrophic currents 

* `Gulf Stream <https://2024-dc-woc-esa.readthedocs.io/en/latest/3_dc_gulfstream/dc_gulfstream_details.html>`_ 
    * Area: [75°W, 45°W, 20°N, 50°N]
    * Dynamical specificities: High variability region with mixed geostrophic and ageostrophic dynamics. 
    * WOC products: Merged SSH/SST currents, Drifter data-driven currents, BFN-QG geostrophic currents, Doppler currents
    
* `Atlantic circumpolar <https://2024-dc-woc-esa.readthedocs.io/en/latest/4_dc_agulhas/dc_agulhas_details.html>`_ 
    * Area: [14°E, 35°E, 45°S, 30°S]
    * Dynamical specificities: Strongly geostrophic region. 
    * WOC products: Drifter data-driven currents, BFN-QG geostrophic currents, Doppler currents
  
:raw-html:`<br />` 
    
.. raw:: html

    <embed> 
        <center>
        <div id="image_map"> <map name="map_example"> <area href="https://2024-dc-woc-esa.readthedocs.io/en/latest/4_dc_agulhas/dc_agulhas_details.html" target="_blank" alt="North Sea" shape=poly coords="335,295, 335,270, 370,270, 370,295"> <area href="https://2024-dc-woc-esa.readthedocs.io/en/latest/3_dc_gulfstream/dc_gulfstream_details.html" target="_blank" alt="Gulf Stream" shape=poly coords="160,150, 160,110, 210,110, 210,150"> <area href="https://2024-dc-woc-esa.readthedocs.io/en/latest/2_dc_medsea/dc_medsea_details.html" target="_blank" alt="Mediterranean Sea" shape=poly coords="300,155, 300,115, 345,115, 345,155">  <img src="https://github.com/ocean-data-challenges/2024_DC_WOC-ESA/assets/33433820/4be7d3ae-ea45-4b5c-b049-06ef985024c0" title="Gulf Stream" alt="image map example" width=600 height=350 usemap="#map_example"></map> </div> </center>

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
   :caption: DC - Mediterranean Sea

   2_dc_medsea/dc_medsea_details.md
   2_dc_medsea/dc_medsea_overalleval.md
   2_dc_medsea/dc_medsea_sshsstproduct.md 
   2_dc_medsea/dc_medsea_drifterproduct.md 
   2_dc_medsea/dc_medsea_bfnqgproduct.md 
   2_dc_medsea/dc_medsea_dopplerproduct.md 
   2_dc_medsea/dc_medsea_otherproducts.md 

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: DC - Gulf Stream
 
   3_dc_gulfstream/dc_gulfstream_details.md
   3_dc_gulfstream/dc_gulfstream_overalleval.md 
   3_dc_gulfstream/dc_gulfstream_sshsstproduct.md 
   3_dc_gulfstream/dc_gulfstream_drifterproduct.md 
   3_dc_gulfstream/dc_gulfstream_bfnqgproduct.md 
   3_dc_gulfstream/dc_gulfstream_dopplerproduct.md 
   3_dc_gulfstream/dc_gulfstream_otherproducts.md 
 

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: DC - Agulhas Current
 
   4_dc_agulhas/dc_agulhas_details.md
   4_dc_agulhas/dc_agulhas_overalleval.md
   4_dc_agulhas/dc_agulhas_drifterproduct.md
   4_dc_agulhas/dc_agulhas_bfnqgproduct.md
   4_dc_agulhas/dc_agulhas_dopplerproduct.md
   4_dc_agulhas/dc_agulhas_otherproducts.md
   
.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Metrics details

   5_metrics_det/metrics_1-standard-insitu.md
   5_metrics_det/metrics_2-effective-resolution.md
   5_metrics_det/metrics_3-position-structures.md
   5_metrics_det/metrics_4-dynamic-characterization.md
   5_metrics_det/metrics_5-lagrangian-diagnostics.md
    

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Scripts

   6_scripts/modules.rst

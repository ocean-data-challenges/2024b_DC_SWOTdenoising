# SSH - Along track metrics

<br>
 

<br>

The ocean surface topography reconstruction is compared with independant data from Saral/AltiKa altimeter. The reconstructed maps are first interpolated onto the independant nadir tracks. The following diagnostics are then performed along these tracks and aggregated in 1° longitude x 1° latitude boxes.

## Statistics


[**Check alongtrack statistics example notebook**](gallery/alongtrack_statistics_description.ipynb)


<br>

### &#x2022; Error variance maps
 

<div style="text-align: justify">
The SSH error variance are computed alongtrack and then averaged in each 1° x 1° box. The resulting error variances are then displayed as maps. Also, a band pass filter can be applied before computing the statistics in order to select specific spatial dynamical scales between 65 km and 500 km. 
</div> 

<br>

![DUACS SSH error variance](figures/Maps_DUACS_errvar_glob.png)
<center> 
  <i>Figure 1: Global DUACS SSH error variance with respect to independent nadir SSH Altika at all scales (left) and scales between 65-500 km (right). </i> 
</center>


<br>

###  &#x2022; Explained variance maps
 
<div style="text-align: justify">
Similarly, the SSH explained variance are computed alongtrack and then averaged in each 1° x 1° box. The resulting explained variances are then displayed as maps. 
</div>

<br>

![DUACS SSH explained variance](figures/Maps_DUACS_explvar_glob.png)
<center> 
  <i>Figure 2: Global DUACS SSH explained variance with respect to independent nadir SSH Altika at all scales (left) and scales between 65-500 km (right). </i> 
</center>

<br>

###  &#x2022; Statistics by regime

<div style="text-align: justify">
The statistics displayed in the previous maps are also averaged in specific geographical areas and regimes (<i>coastal, offshore-high variability, offshore-low variability, equatorial band, Arctic, Antarctic</i>) in order to provide a series of scores summarized in the Leaderboards.
</div>

<br>

| <font size="2">Ocean regime  </font>         |<font size="2"> Methods</font>   |   <font size="2">Err variance score (All scales)</font>  |  <font size="2"> Err variance score (65-500km)</font>  | 
|:----------------------|----------|:---------------------------------:|:-------------------------------:| 
|  <font size="2"> <em>Coastal</em> </font>         |<font size="2"> DUACS </font>    |       <font size="2"> 0.683678  </font>                  |    <font size="2">     0.682788     </font>            |  
| | | | |  
|<font size="2"> <em>Offshore high var</em> </font>  |<font size="2"> DUACS </font>    |      <font size="2">  0.941316 </font>                   |   <font size="2">      0.941316   </font>              |  
| | | | |  
|<font size="2"> <em>Offshore low var</em> </font>   |<font size="2"> DUACS </font>    |      <font size="2">  0.800111 </font>                   |  <font size="2">       0.868348   </font>              |  
| | | | |  
|<font size="2"> <em>Equatorial band</em> </font>   |<font size="2"> DUACS </font>    |     <font size="2">   0.764560  </font>                  |   <font size="2">      0.415881   </font>              |  
| | | | |  
|<font size="2"> <em>Arctic</em>       </font>      |<font size="2"> DUACS </font>    |      <font size="2">  0.667031 </font>               |   <font size="2">      0.585653    </font>             |   
| | | | |  
|<font size="2"> <em>Antarctic</em>   </font>       |<font size="2"> DUACS </font>    |       <font size="2"> 0.415025   </font>                 |   <font size="2">      0.077469    </font>             |  

    
<br>

<br>

## Spectral


[**Check alongtrack spectral example notebook**](gallery/alongtrack_spectral_description.ipynb)

<br> 

### &#x2022; Effective resolution maps

<div style="text-align: justify">
The effective resolution corresponds to the spatiotemporal scales of the features that can be properly resolved in the maps (Ballarotta et al., 2019). The spatiotemporal resolution of the previous level 4 global SLA products was estimated by Chelton et al. (2011, 2014) and Chelton and Schlax (2003) based on estimates of the mapping errors in sea surface height (SSH) fields constructed from altimeter data or spectral ratio analysis between maps and along-track altimeter data. 

To estimate the effective resolution, we first compute the spectral noise-to-signal ratio, i.e. the ratio between the spectral content of the mapping error and the spectral content of independent along-track observations:
</div>

$$\text{NSR}(\lambda_s) = \frac{S_{diff}(\lambda_s)}{S_{obs}(\lambda_s)},$$

where $\lambda_s$ is the spatial wavelength, $S_{diff}(\lambda_s)$ is the power spectral density of the difference (SLA$_{obs}$–SLA$_{map}$) and $S_{obs}(\lambda_s)$ is the spectral density of the independent observation.

<div style="text-align: justify">
Here, we define the effective resolution as the wavelenght at which the noise-to-signal ratio is equal or smaller than 0.5. 
</div>

![DUACS SSH effective resolution](figures/Maps_DUACS_effres_glob.png)
<center> 
  <i>Figure 3: Global DUACS SSH effective resolutions with respect to independent nadir SSH Altika. </i> 
</center>

<br>
 
<div style="text-align: justify">
The maps (see DUACS example in Figure 3) are obtained by averaging the effective resolutions on all available along-track segments in each 1°x1° longitude x latitude box.
</div>

<br>
 
 
## References

- Ballarotta, M., Ubelmann, C., Pujol, M. I., Taburet, G., Fournier, F., Legeais, J. F., ... & Picot, N. (2019). On the resolutions of ocean altimetry maps. Ocean Science, 15(4), 1091-1109.

- Chelton, D. B., & Schlax, M. G. (2003). The accuracies of smoothed sea surface height fields constructed from tandem satellite altimeter datasets. Journal of Atmospheric and Oceanic Technology, 20(9), 1276-1302.

- Chelton, D. B., Schlax, M. G., & Samelson, R. M. (2011). Global observations of nonlinear mesoscale eddies. Progress in oceanography, 91(2), 167-216.

- Chelton, D., Dibarboure, G., Pujol, M. I., Taburet, G., & Schlax, M. G. (2014). The Spatial Resolution of AVISO Gridded Sea Surface Height Fields, OSTST Lake Constance, Germany, 28–31 October 2014.


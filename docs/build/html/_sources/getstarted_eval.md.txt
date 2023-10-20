# Evaluate your own maps

<br> 
 

<br> 

Once you have [installed the data challenge](getstarted_install.md) and [downloaded the data](getstarted_data.md), you can now evaluate your maps. Although, you might want to check first that your sea level anomaly or sea surface currents maps respect a certain format. You can then scroll through the different metrics, checking the DUACS evaluation notebooks as example. 


<br> 

## Maps format

The input netcdf files must contain: 
- `latitude`, `longitude` and `time` dimensions. 
- for sea level anomaly evaluation, a sea level anomaly variable `sla`.
- for currents evaluation, a meridional and a zonal currents variable `ugos` and `vgos`.

Note that some formatting can be done using xarray within the jupyter notebooks. For instance, you can change a variable's name simply with `ds = ds.rename_var({'my_sla':'sla'})`. 

### Sea Level Anomaly maps format

![SLA Format](figures/Maps_format_SLA.png)  

### Currents maps format

![Currents_Format](figures/Maps_format_Currents.png)  

 
<br> 

## Available metrics and where to find them

### SSH evaluation with independant nadir SSH data 
 

#### [Check example 1](https://github.com/ocean-data-challenges/2023a_SSH_mapping_OSE/blob/main/nb_diags_global/ssh_scores_DUACS_geos.ipynb)
 

- **Grid boxes statistics (maps)** 
     
   - **SSH Error variance**      
 
   - **SSH Explained variance** 

 
<br> 

- **Statistics by regimes (scalar scores)**  
 
 
<br> 

- **Spectral effective resolution (maps)**
  


<br> 


<br> 

### Currents evaluation with independant drifter currents 

#### [Check example 2](https://github.com/ocean-data-challenges/2023a_SSH_mapping_OSE/blob/main/nb_diags_global/uv_scores_DUACS_geos.ipynb)
 
- **Grid boxes statistics (maps)**  
     
   - **Currents Error variance**      
 
   - **Currents Explained variance** 

<br> 

- **Statistics by regimes (scalar scores)** 
    

<br> 

- **Zonaly averaged rotary spectra (omega-latitude plots)**  
     

<br> 


<br> 

### Currents evaluation with independant drifter trajectories 

#### [Check example 3](https://github.com/ocean-data-challenges/2023a_SSH_mapping_OSE/blob/main/nb_diags_global/uv_scores_DUACS_geos.ipynb) 

- **Drifter deviation maps at fixed horizons (maps)**

- **Drifter deviation global (horizon series)**

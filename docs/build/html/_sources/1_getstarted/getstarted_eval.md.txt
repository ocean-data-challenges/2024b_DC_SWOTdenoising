# Evaluate your denoising technique

<br> 
 

<br> 

Once you have [installed the data challenge](getstarted_install.md) and [downloaded the data](getstarted_data.md), you can now evaluate your maps. Although, you might want to check first that your sea level anomaly or sea surface currents maps respect a certain format. You can then scroll through the different metrics, checking the DUACS evaluation notebooks as example. 


<br> 

## Data format

The input netcdf files must contain: 
- `num_lines` and `num_pixels` dimensions. 
- `latitude`, `longitude` and `time` variables. 
- for denoising evaluation, a denoised ssh variable `ssha_method_denoised`. 

Note that some formatting can be done using xarray within the jupyter notebooks. For instance, you can change a variable's name simply with `ds = ds.rename_var({'my_sla':'sla'})`. 
  

### Currents maps format
 

 
<br> 

## Available metrics and where to find them
 

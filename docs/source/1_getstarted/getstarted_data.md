# Download the data

<br> 

<br>  

The data are hosted and can be accessed on the MEOM server opendap [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/2024b_DC_SWOTdenoising/catalog.html). You can store the downloaded dataset in the swot_data/ directory. 

To download the SWOT dataset locally, you will need approximately 4 GB of disk space.  

## Data information

The dataset is presented with the following directory structure:

### 1) Data for 1 day repeat cycle 

**swot_1j_share/SWOT_L3_LR_SSH_Expert_500_*_v0.2.nc**
 

### 1) Data for 1 day repeat cycle 

**swot_21j_share/SWOT_L3_LR_SSH_Expert_003_*_v0.1.nc**
 

## Download and read the data

 
``` 
wget https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/fileServer/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/2024b_DC_SWOTdenoising/swot_data.tar.gz 
tar -xvf swot_data.tar.gz  
rm -f swot_data.tar.gz
```

and 


``` 
wget https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/fileServer/meomopendap/extract/MEOM/OCEAN_DATA_CHALLENGES/2024b_DC_SWOTdenoising/denoised_data.tar.gz 
tar -xvf denoised_data.tar.gz  
rm -f denoised_data.tar.gz
```

**Example notebooks**


A notebook to illustrate how to download and read the data is available: [demo_download_swotdata.ipynb](../gallery/read_and_download_regional_data.ipynb)
 
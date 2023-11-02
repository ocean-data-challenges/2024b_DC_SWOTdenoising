import numpy as np 
import xarray as xr 




def retrieve_segments(ds, name_var, len_seg):
    
    var = np.array(ds[name_var])
    var = var[:,~np.isnan(var).all(axis=0)]
    var_ravel = np.ravel(var,'F')
    
    var_seg = list()
    var_segs = list()
    
    for i in range(np.size(var_ravel)):
        
        if ~np.isnan(var_ravel[i]):
            var_seg.append(var_ravel[i])
            if np.size(var_seg) == len_seg: 
                var_segs.append(var_seg)
                var_seg = list()  
        else: 
            var_seg = list()
    
    print('Number of segs',np.shape(var_segs)[0],'of size',np.shape(var_segs)[1])
    
    return var_segs
class RegionInfo:

    def __init__(self,name):
        if name == 'GS_info':
            self.region_info = GS_info
        if name == 'Med_info':
            self.region_info = Med_info
        if name == 'ACC_info':
            self.region_info = ACC_info
        if name == 'Cali_info':
            self.region_info = Cali_info
        if name == 'Mada_info':
            self.region_info = Mada_info
      
    
    

GS_info = {
    'name' : 'GS',                                # region name
    'lon_min' : -55,                            # domain min longitude
    'lon_max' : -45,                            # domain max longitude
    'lat_min' : 30,                             # domain min latitude
    'lat_max' : 40,                              # domain max latitude
    'lon_ticks' : ['55$^\circ$ W','53$^\circ$ W','51$^\circ$ W','49$^\circ$ W','47$^\circ$ W'], # ticks label for plots
    'lat_ticks' : ['30$^\circ$ N','32$^\circ$ N','34$^\circ$ N','36$^\circ$ N','38$^\circ$ N']  # ticks label for plots
}

Med_info = {
    'name' : 'Med',                               # region name
    'lon_min' : -2,                              # domain min longitude
    'lon_max' : 8,                             # domain max longitude
    'lat_min' : 36,                             # domain min latitude
    'lat_max' : 44,                              # domain max latitude
    'lon_ticks' : ['2$^\circ$ W','0$^\circ$','2$^\circ$ E','4$^\circ$ E','6$^\circ$ E'],  # ticks label for plots
    'lat_ticks' : ['36$^\circ$ N','38$^\circ$ N','40$^\circ$ N','42$^\circ$ N'] # ticks label for plots
}


ACC_info = {
    'name' : 'ACC',                              # region name
    'lon_min' : 152,                             # domain min longitude
    'lon_max' : 162,                             # domain max longitude
    'lat_min' : -64,                             # domain min latitude
    'lat_max' : -54,                              # domain max latitude
    'lon_ticks' : ['152$^\circ$ E','154$^\circ$ E','156$^\circ$ E','158$^\circ$ E','160$^\circ$ E'],  # ticks label for plots
    'lat_ticks' : ['64$^\circ$ S','62$^\circ$ S','60$^\circ$ S','58$^\circ$ S','56$^\circ$ S'] # ticks label for plots
}

Mada_info = {
    'name' : 'Mada',                               # region name
    'lon_min' : 44,                              # domain min longitude
    'lon_max' : 54,                             # domain max longitude
    'lat_min' : -12,                             # domain min latitude
    'lat_max' : -2,                              # domain max latitude
    'lon_ticks' : ['44$^\circ$ E','46$^\circ$ E','48$^\circ$ E','50$^\circ$ E','52$^\circ$ E'],  # ticks label for plots
    'lat_ticks' : ['12$^\circ$ S','10$^\circ$ S','8$^\circ$ S','6$^\circ$ S','4$^\circ$ S'] # ticks label for plots
}


Cali_info = {
    'name' : 'Cali',                              # region name
    'lon_min' : -130,                             # domain min longitude
    'lon_max' : -120,                             # domain max longitude
    'lat_min' : 30,                             # domain min latitude
    'lat_max' : 40,                              # domain max latitude
    'lon_ticks' : ['130$^\circ$ W','128$^\circ$ W','126$^\circ$ W','124$^\circ$ W','122$^\circ$ W'],  # ticks label for plots
    'lat_ticks' : ['30$^\circ$ N','32$^\circ$ N','34$^\circ$ N','36$^\circ$ N','38$^\circ$ N'] # ticks label for plots
}
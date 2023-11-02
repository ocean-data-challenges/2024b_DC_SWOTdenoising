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
    'lon_ticks' : ['55$^\circ$ W','53$^\circ$ W','51$^\circ$ W','49$^\circ$ W','47$^\circ$ W'], 
    'lat_ticks' : ['30$^\circ$ N','32$^\circ$ N','34$^\circ$ N','36$^\circ$ N','38$^\circ$ N']
}

Med_info = {
    'name' : 'Med',                               # region name
    'lon_min' : -2,                              # domain min longitude
    'lon_max' : 8,                             # domain max longitude
    'lat_min' : 36,                             # domain min latitude
    'lat_max' : 44                              # domain max latitude
}


ACC_info = {
    'name' : 'ACC',                              # region name
    'lon_min' : 152.,                             # domain min longitude
    'lon_max' : 162.,                             # domain max longitude
    'lat_min' : -64.,                             # domain min latitude
    'lat_max' : -54.                              # domain max latitude
}

Mada_info = {
    'name' : 'Mada',                               # region name
    'lon_min' : 44,                              # domain min longitude
    'lon_max' : 54,                             # domain max longitude
    'lat_min' : -12,                             # domain min latitude
    'lat_max' : -2                              # domain max latitude
}


Cali_info = {
    'name' : 'Cali',                              # region name
    'lon_min' : -130,                             # domain min longitude
    'lon_max' : -120,                             # domain max longitude
    'lat_min' : 30,                             # domain min latitude
    'lat_max' : 40                              # domain max latitude
}
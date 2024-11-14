"""
    **Custom Added Constants**
    
     :r_earth:        Earth radius [m]
     :m_earth:        Earth mass [Kg]
     :c:              speed of light [m/s]
     :Om_T:           Earth rotation speed [rad/s]
     :g:              Gravity constant [m/s^2]
     :SWOT_GND_SPEED: Swot ground speed [m/s]
    
"""

import numpy as np
from scipy.constants import *

r_earth = 6378.137e3      # Earth equat. radius [m]
m_earth = 5.9742e24       # Earth mass [Kg]
c = 2.99792458e8          # speed of light [m/s]
Om_t = 7.2921*1e-5       # Earth rotation speed [rad/s]
g = 9.81                  # Gravity constant [m/s^2]
SWOT_GND_SPEED = 7000     # SWOT ground speed [m/s]

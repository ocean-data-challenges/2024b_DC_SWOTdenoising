from enum import Enum
from typing import List

class ROI(Enum):
    """Define a region of interest.
    
    They are defined with [lon_min, lat_min, lon_max, lat_max]
    in degrees. The ROI also register the LLC4320 faces that cover
    the zone of interest. That have been computes with
    :func:swot.diagnostics_toolbox#faces_covering_area

    See Also
    --------
    :func: swot.diagnostics_toolbox#faces_covering_area
    """
    AGULHAS = ([5.0, -45.0, 45.0, -30.0], [1])
    GLOBAL = ([-180, -90, 180, 90], range(13))
    NEW_CALEDONIA = ([160, -30, 175, -16], [8])

    @property
    def area(self) -> List[float]:
        return self.value[0]
    
    @property
    def llc_faces(self) -> List[int]:
        return self.value[1]
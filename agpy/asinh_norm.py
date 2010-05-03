from matplotlib.colors import Normalize
from matplotlib.cm import cbook
from numpy import ma
import numpy as np

class AsinhNorm(Normalize):
    def __call__(self,value,clip=None):
        if clip is None:
            clip = self.clip

        if cbook.iterable(value):
            vtype = 'array'
            val = ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = ma.array([value]).astype(np.float)

        self.autoscale_None(val)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin==vmax:
            return 0.0 * val
        else:
            if clip:
                mask = ma.getmask(val)
                val = ma.array(np.clip(val.filled(vmax), vmin, vmax),
                                mask=mask)
            result = (ma.arcsinh(val)-np.arcsinh(vmin))/(np.arcsinh(vmax)-np.arcsinh(vmin))
        if vtype == 'scalar':
            result = result[0]
        return result



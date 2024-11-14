import numpy as np
import dask.array as da
import functools
import dataclasses as dc
from typing import Optional
from scipy import signal, ndimage


@dc.dataclass
class Kernel:
    """
    Kernel for convolution filters.

    Parameters
    ----------
    hws
        Half window size
    """

    hws: int

    def compute(self, x):
        """
        Generation function.

        Must be implemented by all subclasses.

        Parameters
        ----------
        x
            Input abscissas of the function
        
        Returns
        -------
        :
            The corresponding kernel values
        """
        pass

    def sample(self):
        """
        Sample a kernel between [-hws, hws]

        Returns
        -------
        :
            Sampled kernel
        """
        x = np.arange(-self.hws, self.hws + 1)
        return self.compute(x)

@dc.dataclass
class Kernel2D:
    """
    Kernel2D composed of two 1D kernels

    Parameters
    ----------
    kernel_x
        Kernel applied along the x dimension
    kernel_y
        Kernel applied along the y dimension
    """
    kernel_x: Kernel
    kernel_y: Optional[Kernel] = None

    def __post_init__(self):
        # Kernel has the same definition along the two dimensions
        if self.kernel_y is None:
            self.kernel_y = self.kernel_x

    def sample(self):
        """
        Sample a 2D kernel.

        Returns
        -------
        :
            The sampled kernel
        """
        kernel_x = self.kernel_x.sample()
        kernel_y = self.kernel_y.sample()
        kernel = np.dot(
            np.expand_dims(kernel_x, 1),
            np.expand_dims(kernel_y, 1).T
        )
        return kernel


@dc.dataclass
class GaussianKernel(Kernel):
    """
    Gaussian kernel.

    Parameters
    ----------
    hws
        Half window size
    sigma
        Standard deviation of the gaussian function
    """
    sigma: float
    
    def compute(self, x):
        y = np.exp(-0.5 * np.square(x / self.sigma))
        return y / y.sum()

@dc.dataclass
class LanczosKernel(Kernel):
    """Lanczos 1D kernel.
    
    Parameters
    ----------
    hws
        Half window size in indexes
    cutoff_frequency
        Cutoff frequency in inverse abscissas steps
    """
    cutoff_frequency: float

    def compute(self, x):
        sinc_filter = 2 * self.cutoff_frequency * np.sinc(2 * self.cutoff_frequency * x)
        sigma_factor = np.sinc(x / self.hws)
        return sinc_filter * sigma_factor


@dc.dataclass
class MeanKernel(Kernel):
    """Mean 1D kernel.
    
    Parameters
    ----------
    hws
        Half window size in indexes
    """

    def compute(self, x):
        return np.ones(x.shape, dtype=float)
    

def convolve_2d_nans(z, /, kernel, mode, cval):
    """Convolution using scipy, ignore invalid values.
    
    Parameters
    ----------
    z
        Input vector of size (T, M, N). 
        The 2 last dimensions are the spatial dimensions
    kernel
        Input kernel of size (1, X, Y). We rely on scipy.ndimage.convolve for the
        multi-dimensional convolution, so all non-spatial dimensions should be 
        set to 1
    
    Returns
    -------
    :
        Filtered array of shape (T, M, N)
    """
    # Padding with 0 before the convolution
    # For each windows, the result will be
    # kernel*x[ii-hws:ii+hws+1] / kernel.sum()
    # This means we do not ignore completely the nan because we should have a weight:
    # kernel[~np.isnan].sum()
    mask = np.isnan(z)
    z[mask] = 0
    z_filtered = ndimage.convolve(z, kernel, mode=mode, cval=cval)

    # Weight correction
    # The following convolution will compute for each window :
    # kernel[~np.isnan].sum() / kernel.sum()
    # Applying this ratio will correct the previous filtering result to completely ignore the nan
    w = np.ones_like(z)
    w[mask] = 0
    correction = ndimage.convolve(w, kernel, mode=mode, cval=cval)

    # Weight correction application and masking
    result = z_filtered / correction
    result[mask] = np.nan
    
    return result

def filter_convolution2d(kernel: Kernel2D, /, z, *, x_axis=1, y_axis=2, mode='reflect', cval=0.0):
    """2D lanczos filtering.
    
    The function is compatible with dask and ignore the invalid values.

    Parameters
    ----------
    kernel
        2D kernel used for the convolution
    z
        Input signal that needs to be filtered. Must be at least 2D
    x_axis
        Position of the X axis
    y_axis
        Position of the Y axis
    mode
        Parameter for scipy.ndimage.convolve for handling the border
    cval
        Parameter for scipy.ndimage.convolve for handling the border
    
    Returns
    -------
    z_filtered: np.ndarray
        Filtered variable with a 2D convolution
    """
    # Transpose the input vector so the core dimensions are at the end of the array
    ndim = z.ndim
    z = np.swapaxes(z, x_axis, ndim - 2)
    if y_axis == ndim - 2:
        z = np.swapaxes(z, x_axis, ndim - 1)
    else:
        z = np.swapaxes(z, y_axis, ndim - 1)

    # 2D kernel computation L(x, y) = L(x)*L(y)
    # Adapt to 3D or 4D if necessary. The X and Y dimensions are at the end
    # of the array
    sampled_kernel = kernel.sample()
    overlap=(kernel.kernel_x.hws, kernel.kernel_y.hws)
    for ii in range(ndim - 2):
        sampled_kernel = sampled_kernel[np.newaxis]
        overlap = (0, *overlap)

    # Convolution for all chunks.
    if len(z.chunks[-2]) > 0 or len(z.chunks[-1]) > 0:
        # If the input array is chunked along the the X or Y dimensions, we need to use map overlap
        result = z.map_overlap(
            convolve_2d_nans,
            overlap,
            kernel=sampled_kernel, mode=mode, cval=cval, dtype=float)
        
        if np.logical_or(
            np.any(np.array(z.chunks[-2]) < overlap[-2]),
            np.any(np.array(z.chunks[-1]) < overlap[-1])):
            print(
                "Warning: map_overlap is used with 'depth' greater than one chunk of the original variable."
                " Additionnal rechunking is added and can decrease performance")
            result = result.rechunk(z.chunks)
    else:
        result = z.map_blocks(convolve_2d_nans, kernel=sampled_kernel, mode=mode, cval=cval, dtype=float)
    
    # Transpose back to the original shape
    result = np.swapaxes(result, ndim - 2, x_axis)
    if x_axis == ndim - 1:
        result = np.swapaxes(result, ndim - 2, y_axis)
    else:
        result = np.swapaxes(result, ndim - 1, y_axis)

    return result

import numpy as np 
import xarray as xr 

def deriv1and2(h):
    gradx_h = gradx(h)
    grady_h = grady(h)
    # gradx_h = np.ma.array(gradx_h, mask = h_derivs.mask, fill_value = 1e9 )
    # grady_h = np.ma.array(grady_h, mask = h_derivs.mask, fill_value = 1e9 )
    grad_h  = gradx_h**2 + grady_h**2

    lap_h = laplacian(h)
    #lap_h = np.ma.array(lap_h, mask = h_derivs.mask, fill_value = 1e9 )

    gradxlap_h = gradx(lap_h)
    gradylap_h = grady(lap_h)
    # gradxlap_h = np.ma.array(gradxlap_h, mask = h_derivs.mask, fill_value = 1e9 )
    # gradylap_h = np.ma.array(gradylap_h, mask = h_derivs.mask, fill_value = 1e9 )
    gradlap_h =  gradxlap_h**2 + gradylap_h**2
    
    return grad_h, lap_h

def gradx(I): 
    """
    Calculates the gradient in the x-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last row is left as 0s.
    """
    
    m, n = I.shape
    M = np.ma.zeros([m,n])

    M[0:-1,:] = np.ma.subtract(I[1::,:], I[0:-1,:])
    return M


def grady(I): 
    """
    Calculates the gradient in the y-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last column is left as 0s.
    """
    
    m, n = I.shape
    M = np.ma.zeros([m,n])
    M[:,0:-1] =  np.ma.subtract(I[:,1::], I[:,0:-1])
    return M


def div(px, py): 
    """
    Calculates the divergence of a 2D field. 
    For the specific application of image denoising, the calculation follows Chambolle (REF)
    ## BELOW, TO BE CLARIFIED
    The x component of M (Mx) first row is = to the first row of px.
    The x component of M (Mx) last row is = to - the before last row of px. (last one = 0)
    The y component of M (My) first column is = to the first column of py.
    The y component of M (My) last column is = to - the before last column of py. (last one = 0)
    ??#(de sorte que div=-(grad)^*)
    Parameters: two 2D ndarray
    Returns: 2D ndarray
    """
    m, n = px.shape
    M = np.ma.zeros([m,n])
    Mx = np.ma.zeros([m,n])
    My = np.ma.zeros([m,n])
 
    Mx[1:m-1, :] = px[1:m-1, :] - px[0:m-2, :]
    Mx[0, :] = px[0, :]
    Mx[m-1, :] = -px[m-2, :]

    My[:, 1:n-1] = py[:, 1:n-1] - py[:, 0:n-2]
    My[:, 0] = py[:,0]
    My[:, n-1] = -py[:, n-2]
     
    M = Mx + My;
    return M


def laplacian(u):
    """
    Calculates laplacian using the divergence and gradient functions defined in the module.
    Parameter: 2D ndarray
    Returns: 2D ndarray
    """
    Ml = div(gradx(u), grady(u));
    return Ml
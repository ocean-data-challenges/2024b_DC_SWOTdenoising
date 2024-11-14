# Gulf Stream region

<br>



## Snapshots
 

|**SSH**| 
|------------| 
| ![1d GS ssh](../_static/plot_intercomp_GS_1d_2023-04-23_ssh.png)|  

|**Gradients of SSH**|
|------------| 
| ![1d GS ssh](../_static/plot_intercomp_GS_1d_2023-04-23_grad.png)|  

|**Laplacian of SSH** |
|------------| 
| ![1d GS ssh](../_static/plot_intercomp_GS_1d_2023-04-23_lapl.png)| 
 
<br>

## Movies
 
 
**SSH**
  
  
<video controls width="800">
  <source src="../_static/movie_intercomp_Gulfstream_ssh.mp4" type="video/mp4" /> 
  
</video>
 
 
 

**Gradients of SSH** 



<video controls width="800">
  <source src="../_static/movie_intercomp_Gulfstream_grad.mp4" type="video/mp4" /> 
  
</video>
  

**Laplacian of SSH**  


<video controls width="800">
  <source src="../_static/movie_intercomp_Gulfstream_lapl.mp4" type="video/mp4" /> 
  
</video>
 
 
<br>

 
## Physical pdf 



| **SSH** | **Gradients of SSH** | **Laplacian of SSH** |
|----|----|----|
| ![1d GS ssh](../_static/pdf_compare_GS_1d_ssh.png) |![1d GS ssh](../_static/pdf_compare_GS_1d_grad.png) | ![1d GS ssh](../_static/pdf_compare_GS_1d_lapl.png) |

 
 
<br>

## Power Spectrum Density 




| **SSH** | **Gradients of SSH** | **Laplacian of SSH** |
|----|----|----|
| ![1d GS ssh](../_static/psd_compare_GS_1d_ssh.png) |![1d GS ssh](../_static/psd_compare_GS_1d_grad.png) | ![1d GS ssh](../_static/psd_compare_GS_1d_lapl.png) |




## Roberts discontinuities

| **SSH** | **Gradients of SSH** |  
|----|----| 
| ![1d GS ssh](../_static/ssh_compare_GS_1d.png) |![1d GS ssh](../_static/grads_compare_GS_1d.png) | 
| **Roberts discontinuities** | **Masked Roberts discontinuities** |  
| ![1d GS ssh](../_static/roberts_compare_GS_1d.png) | ![1d GS ssh](../_static/maskedroberts_compare_GS_1d.png) | 

**Total discontinuity percentages**


Discontinuity scores based on the Roberts discontinuities. The scores correspond to the percentage of yellow points with respect to the blue points in the "Masked Roberts discontinuities" Figure above. 

- Raw SWOT discontinuities: **29.3 %**
- Unet baseline discontinuities: **1.5 %**
- Gomez_V2 discontinuities: **0.4 %** 
- UnetGomez discontinuities: **0.1 %** 

| **Discontinuities function of SWH** |   
|----| 
| ![1d GS ssh](../_static/discontiSWH_compare_GS_1d.png) |

| **Discontinuities spatial distribution** |  
|----| 
|![1d GS ssh](../_static/spatdisconti_compare_GS_1d.png) |
  
# Circular Curve Estimation from Road Centerline String    

EstCureve.py : Estimate cicurlar from Linestring of road centerline. Input road centerline might incloud lead-in/lead-out straight
                lines which come before and after the actural curve.

# Description:

EstCureve.py is a Python script designed to estimate the parameters of a circular curve from a provided linestring representing a road centerline. The script accounts for the potential presence of lead-in/lead-out straight sections preceding the Point of Curve (PC) and following the Point of Tangent (PT). These straight sections can vary in length and are not considered part of the actual circular curve.

The algorithm implemented in EstCureve.py employs a robust curve fitting technique known as RANSAC (RANdom SAmple Consensus) to identify and eliminate potential outliers in the data. This ensures a more accurate estimation of the circular curve parameters even when the input centerline might contain deviations from a perfect circular arc.

# Key Functionalities:

Estimates the center, radius, and other geometric parameters of the circular curve.
Handles lead-in/lead-out straight sections of arbitrary length.
Employs RANSAC to mitigate the influence of outliers in the centerline data.
Benefits:

Provides a reliable method for extracting circular curve information from road centerline data.
Improves the accuracy of curve fitting in scenarios with potential measurement errors or deviations.


_Author : Phisan Santitamnont_ 
          _phisan.chula@gmail.com ; phisan.s@cdg.co.th )_  
         _June,2024_

![REVERSE_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/Round_Abouts.png)

![REVERSE_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/ReverseCurve.png)

![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c0.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c1.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c2.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c3.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c4.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c5.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c6.png)  
![PLOT_CURVE](https://github.com/phisan-chula/CircularCurve/blob/main/CACHE_PrasertManoonkit/Plot_Curve_c7.png)  

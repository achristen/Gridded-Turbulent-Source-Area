# Gridded-Turbulent-Source-Area

Purpose: Calculate and output gridded turbulent source areas following the Kormann and Meixner (2000) analytical source area model. Individual source areas for one time step can be merged into a cumulative source area for a climatology.

## Code

### Developer

#### fpr_kormann_and_meixner.pro

This routine calculates the flux source area ('footprint') in a gridded version based on the following inputs:

required inputs: 
*   z_0_input    : float. roughness length z0 of surface (in m)
*   z_m_input    : float. effective measurement height of flux system (in m) i.e. zm = (z-d)
*   u_input      : float. measured longitudinal wind velocity component (im m/s)
*   sig_v_input  : float. measured standard deviation of lateral wind velocity (im m/s)
*   L_input      : float. measured Obukhov length (in m)
   
optional inputs:
* x_max_input  : float. maximum distance the model grid should extend upwind of the sensor (default 1000 m)
* y_max_input  : float. maximum distance the model grid should extend lateral away from the centreline (default 500 m). Total domain in  y-direction is 2 x y_max_input (default 1000 m)
* d_input      : float. resultion of the grid-cells in m (default 5 m)

The numerical solve used to find the exponents of the power laws for the wind and eddy diffusivity profiles (Eq. 39 & 40 in Kormann & Meixner) works only for a typical range of input parameters. It is possible that no solution is found for a case.

## References

Kormann, R, and Franz X Meixner. 2001. 'An Analytical Footprint Model for Non-Neutral Stratification.' Boundary-Layer Meteorology 99 (2): 207–24.

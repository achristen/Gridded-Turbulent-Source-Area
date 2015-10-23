;+ 
; name: fpr_kormann_and_meixner
; 
; purpose: 
;   calculates the flux source area ('footprint') based on the analytical 
;   footprint model described in Korman and Meixner, Boundary-Layer 
;   Meteorology 99: 207–224 (2001) for given input parameters.
;   
;   #fx_root()# is used to find the exponents of the power laws for the wind
;   and eddy diffusivity profiles (Eq. 39 & 40 in Kormann & Meixner). This 
;   works only for a typical range of input parameters.
;
; category:
;   flux footprint modelling
;
; calling sequence: 
;   footprint = fpr_kormann_and_meixner(z_0_input=z_0_input, $
;       z_m_input=z_m_input, u_input=u_input, sig_v_input=sig_v_input,$
;       L_input=L_input) 
;
; required inputs: 
;   z_0_input    : float. roughness length z0 of surface (in m)
;   z_m_input    : float. effective measurement height of flux system (in m)
;                  i.e. zm = (z-d)
;   u_input      : float. measured longitudinal wind velocity component (im m/s)
;   sig_v_input  : float. measured standard deviation of lateral wind velocity 
;                  (im m/s)
;   L_input      ; float. measured Obukhov length (in m)
;   
; optional inputs:
;   x_max_input  : float. maximum distance the model grid should extend upwind of
;                  the sensor (default 1000 m)
;   y_max_input  : float. maximum distance the model grid should extend lateral 
;                  away from the centreline (default 500 m). Total domain in
;                  y-direction is 2 x y_max_input (default 1000 m)
;   d_input      : float. resultion of the grid-cells in m (default 5 m)
;
; output: 
;   structure    : with tags
;      PHI             DOUBLE  [nx,ny] Flux footprint or vertical flux per unit point source
;      F               DOUBLE  [nx,ny] 
;      D_Y             DOUBLE  [nx,ny]
;      COORD           FLOAT   [nx,ny,2] Geograhical coordinates whith flux system at [0,0].
;                              [nx,ny,0] are x-coordinates for each point of the grid
;                              [nx,ny,1] are y-coordinates for each point of the grid
;      PARAM_M         FLOAT   Exponent of the wind velocity power law
;      PARAM_N         FLOAT   Exponent of the eddy diffusivity power law
;      PARAM_U         DOUBLE  Constant in power-law profile of the wind velocity
;      PARAM_KAPPA     DOUBLE  Constant in power-law profile of the eddy diffusivity
;      PARAM_USTAR     DOUBLE  Friction velocity
;      PARAM_XSI       DOUBLE  Flux length scale
;   
;   to rotate output into wind and output to file, please use function
;   #fpr_write_ncdf#         
;
; subroutines:
;  in the file below
;  #linspace#
;
; reference:
;   Korman and Meixner, Boundary-Layer Meteorology 99: 207–224 (2001).
; 
; example: 
;   footprint = fpr_kormann_and_meixner(z_0_input = 0.02, z_m_input=2.0, $
;   u_input=3.5, sig_v_input=0.35, L_input=1E6)
;   help, footprint, /stru
;   contour, footprint.phi, footprint.coord[*,0,0], footprint.coord[0,*,1], $
;   nlevels=16, /iso, xtitle='Longitudinal (m)', ytitle='Lateral (m)', /ystyle
;   
; revision history: 
;   2006-04-11 kai* - April 11, 2006 in matlab
;   2010-02-15 ac - translated into IDL and implemented different numerical solver
;- 

function fpr_kormann_and_meixner, $
  z_0_input=z_0_input,$
  z_m_input=z_m_input,$
  u_input=u_input,$
  sig_v_input=sig_v_input,$
  L_input=L_input,$
  x_max_input=x_max_input,$
  y_max_input=y_max_input,$
  d_input=d_input
  
  CATCH, Error_status
 
k = 0.4; von-Karman const.

common km_inputs, z_0, z_m, L, z_1, z_2 ; put inputs in common block for subroutines

if not keyword_set(x_max_input) then x_max=1000. else x_max = double(x_max_input)
if not keyword_set(y_max_input) then y_max=500. else y_max = double(y_max_input)
if not keyword_set(d_input) then d=5. else d = double(d_input)
z_0 = double(z_0_input)
z_m = double(z_m_input)
u = double(u_input)
sig_v = double(sig_v_input)
L = double(L_input)

; Bounds of integration for Eq. 42 to 46, defined on p.218
z_1 = 3.*z_0;
z_2 = (1+k)*z_m;

; Implementation of root finding for Eq. 39 & 40
guess_vector = [0.6,0.5,0.1]
m = fx_root(guess_vector,'fequ_m',itmax=1000)
m = real_part(m)

   ;This statement begins the error handler:
   IF Error_status NE 0 THEN BEGIN
      message, 'Model convergence error: ', /informational
      message, 'Error message: ', /informational
      goto, end_of_program
      CATCH, /CANCEL
   ENDIF

guess_vector = [0.0,0.5,1.5]
n = fx_root(guess_vector,'fequ_n',itmax=1000)
n = real_part(n)

   ;This statement begins the error handler:
   IF Error_status NE 0 THEN BEGIN
      message, 'Model convergence error: ', /informational
      message, 'Error message: ', /informational
      goto, end_of_program
      CATCH, /CANCEL
   ENDIF

z_m = double(z_m_input)

; Inversion of Eq. 31
u_star = u * k / (alog(z_m/z_0) + psi_m(z_m,L));

; Eq (41)
U = u_star/k * (I_2(m,z_0/z_m,z_1,z_2,z_m) + J_1(m,z_1,z_2,z_m,L,'f_psi_m')) $
    / (I_1(2.*m,z_1,z_2,z_m) * z_m^m);

; Eq (41)
kappa = k*u_star * J_1(n,z_1,z_2,z_m,L,'f_z_phi_c') $
    / (I_1(2.*n,z_1,z_2,z_m) * z_m^(n-1.));

; Definition of the grid over which the footprint phi is going to be
; calculated
xs = findgen(x_max/d)*d;[1:d:x_max];
ys = (findgen(2*y_max/d+1) - (y_max/d))*d;[-y_max:d:y_max-1];
x = fltarr(n_elements(xs),n_elements(ys),2)
for j=0, n_elements(ys)-1 do x[*,j,0] = xs[*];
for j=0, n_elements(xs)-1 do x[j,*,1] = ys[*];

; Cross wind integrated footprint
; r is defined at the top of p.213, mu after Eq. 18
r = 2. + m - n
mu  = (1. + m ) / r
; Eq. 19
xsi = U * z_m^r / (r^2 * kappa);
; Eq. 21
f = gamma(mu)^(-1) * xsi^mu/(x[*,*,0]^(1+mu)) * exp(-xsi/x[*,*,0]);

; Cross wind diffusion
; Eq. 18
u_bar = gamma(mu)/gamma(1./r) * (r^2*kappa/U)^(m/r)*U*x[*,*,0]^(m/r);
; Eq. 9, definition of sig right after it
sig = sig_v*x[*,*,0]/u_bar;
D_y = (sqrt(2*!pi)*sig)^(-1) * exp(-x[*,*,1]^2/(2.*sig^2));

; Eq. 8
phi = D_y * f;

; Output in gridded form
phi = reform(phi*d^2,n_elements(xs),n_elements(ys));
D_y = reform(D_y*d,n_elements(xs),n_elements(ys));
f   = reform(f*d,n_elements(xs),n_elements(ys));

model = {phi:phi,f:f,D_y:D_y,coord:x,param_m:m,param_n:n,param_u:u,param_kappa:kappa,param_ustar:u_star,param_xsi:xsi}

return, model

end_of_program:

return, -1

end

function psi_m, z, L;
 z = float(z)
 L = float(L)
 zeta = (1. - 16. * float(z)/L)^(1./4);
 if L gt 0 then begin
    phi_c = 5.*z/L;
 endif else begin
    phi_c = (-2.)*alog((1.+zeta)/2) - alog((1.+zeta^2)/2) + 2. * atan(zeta) - !pi/2;
 endelse
 return, phi_c
end

function phi_c, z, L ;
 z = float(z)
 L = float(L)
 if L gt 0 then begin
    phi_c = 1. + 5.*z/L;
 endif else begin
    phi_c = (1. - (16. * z/L))^(-1./2);
 endelse
 return, phi_c
end

function I_1, p, z_1, z_2, z_m
 ;% Eq (42)
 z = linspace(z_1/z_m,z_2/z_m,1000);
 dz = diff(z);
 z = z[0:n_elements(z)-2] + dz/2;
 d = total( z^p *dz );
 return, d
end 

function I_2, p, z0, z_1, z_2, z_m
 ;% Eq (43)
 z = linspace(z_1/z_m,z_2/z_m,1000);
 dz = diff(z);
 z = z[0:n_elements(z)-2]  + dz/2;
 d = total( z^p * alog(z/z0)*dz );
 return, d
end

function I_3, p, z0, z_1, z_2, z_m
 ;% Eq (44)
 z = linspace(z_1/z_m,z_2/z_m,1000);
 dz = diff(z);
 z = z[0:n_elements(z)-2]  + dz/2;
 d = total( z^p * alog(z) * alog(z/z0)*dz );
 return, d
end

function J_1, p, z_1, z_2, z_m, L, f
 ;% Eq (45)
 z = linspace(z_1/z_m,z_2/z_m,1000);
 dz = diff(z);
 z = z[0:n_elements(z)-2]  + dz/2;
 d = total( z^p * CALL_FUNCTION(f,z*z_m,z_m,L) *dz);
 return, d
end

function J_2, p, z_1, z_2, z_m, L, f
 ;% Eq (46)
 z = linspace(z_1/z_m,z_2/z_m,1000);
 dz = diff(z);
 z = z[0:n_elements(z)-2]  + dz/2;
 d = total( z^p * CALL_FUNCTION(f,z*z_m,z_m,L) *alog(z) * dz );
 return, d
end

;% Functions that are the arguments of J_1 and J_2 

function f_psi_m, z, z_m, L
 z = float(z)
 d = psi_m(z,L)
 return, d
end

function f_z_phi_c, z, z_m, L
 z = float(z)
 d = z/(phi_c(z,L)*z_m)
 return, d
end

function fequ_m, m
 common km_inputs, z_0, z_m, L, z_1, z_2
 ;% Eq (39)
 all_d = fltarr(n_elements(m))
 for j=0, n_elements(m)-1 do begin
  A = I_1(2.*m[j],z_1,z_2,z_m)   * ( I_3(m[j],z_0/z_m,z_1,z_2,z_m) + J_2(m[j],z_1,z_2,z_m,L,'f_psi_m') );
  B = I_2(2.*m[j],1,z_1,z_2,z_m) * ( I_2(m[j],z_0/z_m,z_1,z_2,z_m) + J_1(m[j],z_1,z_2,z_m,L,'f_psi_m') );
  all_d[j] = B - A;
 endfor
 return, all_d
end

function fequ_n, n
 common km_inputs, z_0, z_m, L, z_1, z_2
 all_d = fltarr(n_elements(n))
 for j=0, n_elements(n)-1 do begin
  ;% Eq (40)
  A = I_1(2.*n[j],z_1,z_2,z_m)   * J_2(n[j],z_1,z_2,z_m,L,'f_z_phi_c')
  B = I_2(2.*n[j],1,z_1,z_2,z_m) * J_1(n[j],z_1,z_2,z_m,L,'f_z_phi_c')
  all_d[j] = B - A
 endfor
 return, all_d
end

function diff, array
  return, array[1:n_elements(array)-1] - array[0:n_elements(array)-2]
end

;; THE FOLLOWING ROUTINES HAVE BEEN CODED EXTRERNALLY

;;;;;;;;;;;;;;;;;;;;;;;;
;;;   linspace.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   21-Jun-2001
;;;  Version: 0.35 (CVS: $Revision: 1.3 $)
;;;  Description:
;;;     Return a real vector of length N, containing equidistant values
;;;   between x1 and x2 inclusively.
;;;   If N is omitted, a value of 100 is adopted.
;;;     Mimics the octave function `linspace'.
;;;     Allows for the first 2 arguments to be written as a vector,
;;;   so `linspace(minmax(v),10)' works.
;;;  Keywords:
;;;   PERIODIC  -- flag for periodic grid (i.e. x[n-1]=x2-dx)
;;;   GHOST     -- set this to the number of ghost cells before x1 or
;;;                after x2; GHOST can be a 2-element array
;;;                [left,right] or a scalar (applied to both sides)
;;;   UNIQUE    -- flag for returning a list of unique elements.
;;;                This implies that you may get less than N elements
;;;                (in many cases just one).
;;;                Useful if you call
;;;                  contour,..,LEVELS=linspace(minmax(var),N)
;;;                where no repeated levels are allowed

function linspace, x1, x2, n, $
  PERIODIC=peri, $
  GHOST=ghost, $
  UNIQUE=unique
  on_error,2

  if (n_elements(ghost) eq 0) then begin
    nghost = [0,0]
  endif else begin
    nghost = rebin([ghost], 2)  ; replicate if necessary
  endelse
  if (keyword_set(peri)) then nghost=[0,-1]
  ;default, unique, 0

  if (n_elements(x1) ge 2) then begin ; Reshuffle arguments
    xx1 = x1[0]
    xx2 = x1[1]
    if (n_elements(x2) ne 0) then nn = x2
  endif else begin
    xx1 = x1
    xx2 = x2
    if (n_elements(n) ne 0) then nn = n
  endelse
  ;default, nn, 100
  n_real = nn - nghost[0] - nghost[1]
  list = xx1 + (findgen(nn)-nghost[0])*(xx2-xx1)/(n_real-1)

  if keyword_set(unique) then list = list(uniq(list))

  return, list
end

; End of file linspace.pro
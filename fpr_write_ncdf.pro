;+ 
; name: fpr_write_ncdf
; 
; purpose: 
;   calculates the flux source area ('footprint') based on available 
;   footprint models (currently only Korman and Meixner), rotates output 
;   into mean wind and writes a geographcally referenced raster into
;   a NCDF file.
;   
;   #fpr_korman_and_meixner# is used to calculate the footprints
;
; category:
;   flux footprint modelling
;
; calling sequence: 
;   fpr_write_ncdf, filename = filename, z_0_input=z_0_input, $
;    z_m_input=z_m_input, u_input=u_input, wd_input=wd_input, $
;    sigv_input=sigv_input, L_input=L_input, domain_x=domain_x, $
;    domain_y=domain_y, domain_output=domain_output, $
;    d_input=d_input, juliantime = juliantime, site=site, $
;    timezone=timezone, provider=provider
;
; required inputs: 
;   z_0_input    : float. roughness length z0 of surface (in m)
;   z_m_input    : float. effective measurement height of flux system (in m)
;                  i.e. zm = (z-d)
;   u_input      : float. measured longitudinal wind velocity component (im m/s)
;   sig_v_input  : float. measured standard deviation of lateral wind velocity 
;                  (im m/s)
;   L_input      : float. measured Obukhov length (in m)
;   juliantime   : double. time of the footprint
;   wd_input     ; float. wind directions in degree from geographic North.
;   
; optional inputs:
;   filename     : string. path. filename of the ncdf file the footprint will 
;                  be written to.
;   domain_output; float. the domain size in m for the ncdf file, where the
;                  flux system will be in the center (i.e. domain size will be
;                  domain_output x domain_output)
;   x_max_input  : float. maximum distance the model grid should extend upwind of
;                  the sensor (default 1000 m)
;   y_max_input  : float. maximum distance the model grid should extend lateral 
;                  away from the centreline (default 500 m). Total domain in
;                  y-direction is 2 x y_max_input (default 1000 m)
;   d_input      : float. resultion of the grid-cells in m (default 5 m)
;   site         : string. name of site / system
;   timezone     : string. time zone of time information.
;   provider     : data providor or operator of site.
;   
; output: 
;   NCDF file with all inputs and outputs of the model.     
;
; subroutines:
;  #fpr_korman_and_meixner#
;
; example: 
;  fpr_write_ncdf, filename='/Users/achristn/Desktop/footprint.nc', z_0_input=0.2, $
;  z_m_input=10.0, u_input=3.0, wd_input=45.0, sigv_input=0.5, L_input=1E6, $
;  juliantime=julday(1,1,2000)
;   
; revision history: 
;   2010-02-15 ac
;- 

pro fpr_write_ncdf, $
  $
  filename = filename, $
  z_0_input=z_0_input, $
  z_m_input=z_m_input, $
  u_input=u_input, $
  wd_input=wd_input, $
  sigv_input=sigv_input, $
  L_input=L_input, $
  domain_x=domain_x, $
  domain_y=domain_y, $
  domain_output=domain_output, $
  d_input=d_input, $
  juliantime = juliantime, $
  site=site, $
  timezone=timezone, $
  provider=provider
  
  ;*****************************
  ; handle command line input
  ;*****************************
  
  command_line_input = COMMAND_LINE_ARGS(COUNT=commandline_count)
  
  if commandline_count gt 0 then begin
    
    if commandline_count gt 6 then begin
    filename = command_line_input[0]
    message, 'Filename: '+string(filename), /informational
    z_0_input = float(command_line_input[1])
    message, 'Roughness length (m): '+string(z_0_input), /informational
    z_m_input = float(command_line_input[2])
    message, 'Measurement height (m): '+string(z_m_input), /informational
    u_input = float(command_line_input[3])
    message, 'Mean wind speed (m/s): '+string(u_input), /informational
    wd_input = float(command_line_input[4])
    message, 'Wind direction (degrees): '+string(wd_input), /informational
    sigv_input = float(command_line_input[5])
    message, 'Standard deviation of lateral wind (m/s): '+string(sigv_input), /informational
    L_input = float(command_line_input[6])
    message, 'Okuhov length (m): '+string(L_input), /informational
    
    ; optional inputs
    
    if commandline_count gt 7 then begin
      domain_x = float(command_line_input[7])
      message, 'Domain size for calculation X (m): '+string(domain_x), /informational
    endif
    if commandline_count gt 8 then begin
      domain_y = float(command_line_input[8])
      message, 'Domain size for calculation Y (m): '+string(domain_y), /informational
    endif
    if commandline_count gt 8 then begin
      domain_output = float(command_line_input[9])
      message, 'Domain size for file output (m): '+string(domain_output), /informational
    endif
    if commandline_count gt 10 then begin
      d_input = float(command_line_input[10])
      message, 'Grid resolution (m): '+string(d_input), /informational
    endif
    if commandline_count gt 11 then begin
      juliantime = double(command_line_input[11])
      message, 'Julian time of dataset: '+string(juliantime), /informational
    endif
    if commandline_count gt 12 then begin
      site = string(command_line_input[12])
      message, 'Site name: '+string(site), /informational
    endif
    if commandline_count gt 13 then begin
      site = string(command_line_input[13])
      message, 'Time zone: '+string(site), /informational
    endif
    if commandline_count gt 14 then begin
      provider = string(command_line_input[14])
      message, 'Data Provider: '+string(provider), /informational
    endif
    endif else message, 'Not enough arguments provided. Need at least 6 required arguments supplied.'

  endif
  
  ;*****************************
  ; check regular inputs
  ;*****************************
  
  z_0_input = float(z_0_input[0])
  z_m_input = float(z_m_input[0])
  u_input = float(u_input[0])
  wd_input = float(wd_input[0]) 
  sigv_input = float(sigv_input[0])
  L_input = float(L_input[0])
  
  if not keyword_set(d_input) then d_input = 5.0
  if not keyword_set(domain_output) then domain_output = 2000.0
  if not keyword_set(domain_x) then domain_x = 2000.0
  if not keyword_set(domain_y) then domain_y = 500.0
  if not keyword_set(site) then site = 'not defined'
  if not keyword_set(timezone) then timezone = 'not defined'
  if not keyword_set(provider) then provider = 'not defined'
  title = 'Flux footprint model output'
  
  ;*****************************
  ; run model
  ;*****************************
  
  footprint = fpr_kormann_and_meixner( $
    z_0_input=z_0_input,$
    z_m_input=z_m_input,$
    u_input=u_input,$
    sig_v_input=sigv_input,$
    L_input=L_input,$
    x_max_input=domain_x,$
    y_max_input=domain_y,$
    d_input=d_input)
  
  return_analysis = size(footprint)
  
  if return_analysis[2] eq 8 then begin
  
  ;*****************************
  ; extract distance of max. footprint
  ;*****************************
  
  ixs = where(footprint.f eq max(footprint.f,/NaN))
  x_distance = footprint.coord[*,*,0]
  x_max = x_distance[ixs[0]]
  
  ;*****************************
  ; pad and rotate footprint
  ;*****************************
  
  ;fill it up to a 2 by 2 km raster
  fx = 1+domain_output*2./d_input & fy = fx
  full = dindgen(fx,fy) & full[*] = 0.0 ; double size of output area before rotation
  full[(fx-1)/2:(fx-1)/2+(domain_x/d_input)-1,(fy/2-1)-domain_y/d_input:(fy/2-1)+domain_y/d_input] = footprint.phi
  undef = where(finite(full) eq 0, ucnt)
  if ucnt gt 0 then full[undef] = 0.0 ; set to zero
  rotang = wd_input-90.
  rotated = rot(full,rotang,1,(fx-1)/2,(fy-1)/2,/pivot)
  rotated = rotated[(fx-1)/2-domain_output/(2*d_input):(fx-1)/2+domain_output/(2*d_input),(fy-1)/2-domain_output/(2*d_input):(fy-1)/2+domain_output/(2*d_input)] ; cut subset
  coordinates = findgen(domain_x/d_input+1) * d_input - domain_x/2
  nx = n_elements(rotated[*,0])
  ny = n_elements(rotated[0,*])
  
  ;*****************************
  ; write to ncdf file
  ;*****************************
  
  ;set-up file  -----------------------------------
  
  if not keyword_set(filename) then filename=dialog_pickfile(default_extension='nc',file='footprint',/overwrite_prompt,/write,title='Save Footprint as...')
  
  caldat, juliantime, mon, day, year, hour, minute
  strtime=strupcase(string(year,format='(i4.4)')+'-'+string(mon,format='(i2.2)')+'-'+string(day,format='(i2.2)')+' '+string(hour,format='(i2.2)')+':'+string(minute,format='(i2.2)'))
  filewritten=strupcase(systime())
  description = 'Turbulent flux footprint calculated based on Korman and Meixner, Boundary-Layer Meteorology 99: 207–224 (2001) for '+strtime+' ('+timezone+') at '+site+'.'
  
  ;create new file / overwrite existing file ----------------

  cdfid  = NCDF_CREATE(filename, /clobber)

  ;define dimensions ----------------------------------------

  dim_x = NCDF_DIMDEF(cdfid,'EASTING',nx)
  dim_y = NCDF_DIMDEF(cdfid,'NORTHING',ny)
  dim_i = NCDF_DIMDEF(cdfid,'MODEL_INPUT',1)
  dim_p = NCDF_DIMDEF(cdfid,'MODEL_PARAM',1)
   
   ;define global attributes .--------------------------------

  NCDF_ATTPUT, cdfid, /GLOBAL, 'TITLE', title, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'DESCRIPTION', description, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'PROVIDER', provider, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'SITE', site, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'TIME', strtime, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'JULIANTIME', juliantime, /double
  NCDF_ATTPUT, cdfid, /GLOBAL, 'TIMEZONE', timezone, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'FILE_CREATED', filewritten, /char
   
   ;define variables ------------------------------------------

  id_x  = NCDF_VARDEF(cdfid,'EASTING',dim_x,/float)
  id_y  = NCDF_VARDEF(cdfid,'NORTHING',dim_y,/float) 
  id_phi= NCDF_VARDEF(cdfid,'PHI',[dim_x,dim_y],/double) ; footprint
  
  id_d  = NCDF_VARDEF(cdfid,'MODEL_INPUT_D',dim_i,/double)
  id_z0 = NCDF_VARDEF(cdfid,'MODEL_INPUT_Z0',dim_i,/float) ; zero-plane displacement
  id_zm = NCDF_VARDEF(cdfid,'MODEL_INPUT_ZM',dim_i,/float) ; effective measurement height
  id_L  = NCDF_VARDEF(cdfid,'MODEL_INPUT_OBUKHOV_LENGTH',dim_i,/float) ; obukhov length
  id_sv = NCDF_VARDEF(cdfid,'MODEL_INPUT_SIGMA_V',dim_i,/float) ; obukhov length
  id_u  = NCDF_VARDEF(cdfid,'MODEL_INPUT_WIND_VELOCITY',dim_i,/float) ; obukhov length
  id_wd = NCDF_VARDEF(cdfid,'MODEL_INPUT_WIND_DIRECTION',dim_i,/float) ; wind direction
   
  id_pm = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_M',dim_p,/float) ; Exponent of the wind velocity power law
  id_pn = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_N',dim_p,/float) ; Exponent of the eddy diffusivity power law
  id_uc = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_CONST_U',dim_p,/float) ; Constant in power-law profile of the wind velocity
  id_us = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_USTAR',dim_p,/float) ; Friction velocity
  id_ka = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_KAPPA',dim_p,/float) ;Constant in power-law profile of the eddy diffusivity
  id_xi = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_XSI',dim_p,/float) ;Flux length scale
  id_xm = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_XMAX',dim_p,/float) ;Distance of maximum flux per unit point source 
  
   ;define variable attributes --------------------------------
   
  NCDF_ATTPUT, cdfid, id_x,  'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_y,  'UNIT', 'metre', /char  
  NCDF_ATTPUT, cdfid, id_phi,'UNIT', 'metre^-2', /char 
  NCDF_ATTPUT, cdfid, id_phi,'DESCRIPTION', 'Flux footprint or vertical flux per unit point source'
  NCDF_ATTPUT, cdfid, id_d,  'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_d,  'DESCRIPTION', 'Spatial resolution of flux footprint raster data'
  NCDF_ATTPUT, cdfid, id_z0, 'UNIT', 'metre', /char
  NCDF_ATTPUT, cdfid, id_z0, 'DESCRIPTION', 'Roughness length', /char 
  NCDF_ATTPUT, cdfid, id_zm, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_zm, 'DESCRIPTION', 'Effective measurement height (z-d)', /char 
  NCDF_ATTPUT, cdfid, id_L , 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_L,  'DESCRIPTION', 'Obukhov length at zm', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Standard deviation of lateral wind component at zm', /char
  NCDF_ATTPUT, cdfid, id_u , 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Longitudinal wind velocity component at zm', /char
  NCDF_ATTPUT, cdfid, id_wd, 'UNIT', 'degree', /char
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Wind direction relative to geographic North (clockwise) at zm', /char
  NCDF_ATTPUT, cdfid, id_pm, 'UNIT', '-', /char  
  NCDF_ATTPUT, cdfid, id_pm, 'DESCRIPTION', 'Exponent of the wind velocity power law', /char  
  NCDF_ATTPUT, cdfid, id_pn, 'UNIT', '-', /char 
  NCDF_ATTPUT, cdfid, id_pn, 'DESCRIPTION', 'Exponent of the eddy diffusivity power law', /char 
  NCDF_ATTPUT, cdfid, id_uc, 'UNIT', 'metre^(1−m) second^−1', /char 
  NCDF_ATTPUT, cdfid, id_uc, 'DESCRIPTION', 'Constant in power-law profile of the wind velocity', /char 
  NCDF_ATTPUT, cdfid, id_us, 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_us, 'DESCRIPTION', 'Friction velocity', /char 
  NCDF_ATTPUT, cdfid, id_ka, 'UNIT', 'metre^(2−n) second^−1', /char 
  NCDF_ATTPUT, cdfid, id_ka, 'DESCRIPTION', 'Constant in power-law profile of the eddy diffusivity', /char 
  NCDF_ATTPUT, cdfid, id_xi, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_xi, 'DESCRIPTION', 'Flux length scale', /char 
  NCDF_ATTPUT, cdfid, id_xm, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_xm, 'DESCRIPTION', 'Upwind distance of maximum flux per unit point source', /char 
  
  ;put the file in data mode ---------------------------------

  NCDF_CONTROL, cdfid, /endef

  ;put data in variable --------------------------------------
  
  NCDF_VARPUT, cdfid, id_x, coordinates
  NCDF_VARPUT, cdfid, id_y, coordinates
  NCDF_VARPUT, cdfid, id_phi, rotated
  NCDF_VARPUT, cdfid, id_d, d_input
  NCDF_VARPUT, cdfid, id_z0, z_0_input
  NCDF_VARPUT, cdfid, id_zm, z_m_input
  NCDF_VARPUT, cdfid, id_L, L_input
  NCDF_VARPUT, cdfid, id_sv, sigv_input
  NCDF_VARPUT, cdfid, id_u, u_input
  NCDF_VARPUT, cdfid, id_wd, wd_input
  
  NCDF_VARPUT, cdfid, id_pm, footprint.param_m
  NCDF_VARPUT, cdfid, id_pn, footprint.param_n
  NCDF_VARPUT, cdfid, id_uc, footprint.param_u
  NCDF_VARPUT, cdfid, id_us, footprint.param_ustar
  NCDF_VARPUT, cdfid, id_ka, footprint.param_kappa
  NCDF_VARPUT, cdfid, id_xi, footprint.param_xsi
  NCDF_VARPUT, cdfid, id_xm, x_max
  
  ;close ncdf file -------------------------------------------

  NCDF_CLOSE, cdfid
  message, 'Successfully written netCDF file of footprint to '+filename, /informational
  
  endif else message, 'Model error - Unable to write file.', /informational

end
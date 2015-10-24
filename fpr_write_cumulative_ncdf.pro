;+
; name: fpr_write_cumulative_ncdf
;
; purpose:
;   aggregates a cumulative flux source area ('footprint') file by averaging
;   multiple individual footprint files form multiple time steps. The individual 
;   files must have been generated in netCDF format by #fpr_write_ncdf# and have 
;   the same dimensions. By selectively choosing files, one can create cumulative
;   footprints for specific cases (e.g. night, day).
;
; category:
;   flux footprint modelling
;
; calling sequence:
;   fpr_write_cumulative_ncdf, ncfiles, outfile=outfile, description=description
;   
;   "ncfiles". stringarray. 
;   This is an optional parameter that supplies of a string array of paths to all 
;   individual foortprint files that should be aggregated. If not provided, a file
;   dialog will promt the user to select multiple footprint files.
;   
;   "outfile". string (path).
;   This is the file path to which the cumulative file (to be created) will be written
;   to. If not provided, a file dialog will promt the user to choose the name and location
;   of the cumulative footprint output.
;   
;   "description". string.
;   This is a used-defined, optional description that allows to place a description into
;   the netCDF output, e.g. "night-time cases".
;
; output:
;   A single NCDF file that contains the cumulative footprint. The file hsas the
;   same structure as individual footprint files in netCDF format.
;
; subroutines:
;  #fpr_write_ncdf#
;  #ncdf2struct# (contained at end of this file)
;
; example:
;
; revision history:
;   2010-02-15 ac
;-

pro fpr_write_cumulative_ncdf, ncfiles, outfile=outfile, description=description

  if n_params() lt 1 then ncfiles=dialog_pickfile(default_extension='nc',/read,title='Select multiple Footprint netCDFs...', /multiple, path='/Users/Shared/Footprint_Model_Data/')
  if not keyword_set(outfile) then filename=dialog_pickfile(default_extension='nc',file='footprint',/overwrite_prompt,/write,title='Save Cummulative footprint as...')
  if not keyword_set(description) then description=''
  
  cummulative_phi = 0D
  nf = n_elements(ncfiles)

  for i=0, nf-1 do begin
    message, 'Reading: '+ncfiles[i], /informational
    footprint = ncdf2struct(ncfiles[i])
    cummulative_phi = cummulative_phi + double(footprint.phi) / nf
    if i eq 0 then begin
       juliantime = footprint.juliantime
       z_0_input = footprint.model_input_z0
       z_m_input = footprint.model_input_zm
       l_input = footprint.model_input_obukhov_length
       sigv_input = footprint.model_input_sigma_v
       u_input = footprint.model_input_wind_velocity
       wd_input = footprint.model_input_wind_direction
       d_input = footprint.model_input_d
       xsi = footprint.model_parameter_xsi
       x_max = footprint.model_parameter_xmax
       param_m = footprint.model_parameter_m
       param_n = footprint.model_parameter_n
       param_u = footprint.model_parameter_const_u
       ustar = footprint.model_parameter_ustar
       param_kappa = footprint.model_parameter_kappa  
       site = footprint.site
       timezone = footprint.timezone
       provider = footprint.provider
       easting = footprint.easting
       northing = footprint.northing
       total_phi = total(footprint.phi)
    endif else begin
       juliantime = [juliantime,footprint.juliantime]    
       z_0_input = [z_0_input,footprint.model_input_z0]
       z_m_input = [z_m_input,footprint.model_input_zm]
       l_input = [l_input,footprint.model_input_obukhov_length]
       sigv_input = [sigv_input,footprint.model_input_sigma_v]
       u_input = [u_input,footprint.model_input_wind_velocity] 
       wd_input = [wd_input,footprint.model_input_wind_direction] 
       x_max = [x_max,footprint.model_parameter_xmax] 
       xsi = [xsi,footprint.model_parameter_xsi]
       param_m = [param_m,footprint.model_parameter_m]
       param_n = [param_n,footprint.model_parameter_n]
       param_u = [param_u,footprint.model_parameter_const_u]
       ustar = [ustar,footprint.model_parameter_ustar]
       param_kappa = [param_kappa,footprint.model_parameter_kappa]
       total_phi = [total_phi,total(footprint.phi)]
    endelse
    if strupcase(site) ne strupcase(footprint.site) then message, 'Error: Input files from different sites.'
    if strupcase(timezone) ne strupcase(footprint.timezone) then message, 'Error: Input files with different time zones.'
    for j = 0, n_elements(easting)-1 do begin
      if easting[j] ne footprint.easting[j] then begin
        message, 'Error: Input files with different grid (easting).'
        stop
      endif 
    endfor 
    for j = 0, n_elements(northing)-1 do begin
      if northing[j] ne footprint.northing[j] then begin
        message, 'Error: Input files with different grid (northing).'
        stop
      endif  
    endfor
  endfor  
  nx = n_elements(footprint.easting)
  ny = n_elements(footprint.northing)
  
  ;**********************************************
  ; write cummulative footprint to ncdf file
  ;**********************************************
  
  ;set-up file  -----------------------------------
  
  caldat, min(juliantime), min_mon, min_day, min_year, min_hour, min_minute
  caldat, max(juliantime), max_mon, max_day, max_year, max_hour, max_minute
  
  min_strtime=strupcase(string(min_year,format='(i4.4)')+'-'+string(min_mon,format='(i2.2)')+'-'+string(min_day,format='(i2.2)')+' '+string(min_hour,format='(i2.2)')+':'+string(min_minute,format='(i2.2)'))
  max_strtime=strupcase(string(max_year,format='(i4.4)')+'-'+string(max_mon,format='(i2.2)')+'-'+string(max_day,format='(i2.2)')+' '+string(max_hour,format='(i2.2)')+':'+string(max_minute,format='(i2.2)'))

  str_min_time=strupcase(min_strtime)
  str_max_time=strupcase(max_strtime)
  filewritten=systime()
  description_gen = 'Cumulative turbulent flux footprint calculated based on Korman and Meixner, Boundary-Layer Meteorology 99: 207–224 (2001) at '+site+'.'
  
  title = 'Cumulative flux footprint'
  
  ;create new file / overwrite existing file ----------------

  cdfid  = NCDF_CREATE(outfile, /clobber)

  ;define dimensions ----------------------------------------

  dim_x  = NCDF_DIMDEF(cdfid,'EASTING',nx)
  dim_y  = NCDF_DIMDEF(cdfid,'NORTHING',ny)
  dim_i  = NCDF_DIMDEF(cdfid,'GLOBAL_MODEL_INPUT',1)
  dim_ii = NCDF_DIMDEF(cdfid,'MODEL_INPUT_RUN_INDIVIDUAL',n_elements(z_0_input))
  dim_p  = NCDF_DIMDEF(cdfid,'GLOBAL_MODEL_PARAM',1)
  dim_pp = NCDF_DIMDEF(cdfid,'MODEL_PARAM_RUN_INDIVIDUAL',n_elements(z_0_input))
  
   ;define global attributes .--------------------------------

  NCDF_ATTPUT, cdfid, /GLOBAL, 'TITLE', title, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'DESCRIPTION', description_gen, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'DESCRIPTION_TIMEPERIOD', description, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'PROVIDER', provider, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'SITE', site, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'TIME_FROM', str_min_time, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'TIME_TO', str_max_time, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'JULIANTIME', juliantime, /double
  NCDF_ATTPUT, cdfid, /GLOBAL, 'TIMEZONE', timezone, /char
  NCDF_ATTPUT, cdfid, /GLOBAL, 'FILE_CREATED', filewritten, /char
   
   ;define variables ------------------------------------------

  id_x  = NCDF_VARDEF(cdfid,'EASTING',dim_x,/float)
  id_y  = NCDF_VARDEF(cdfid,'NORTHING',dim_y,/float) 
  id_phi= NCDF_VARDEF(cdfid,'PHI',[dim_x,dim_y],/double) ; cummulative footprint
  
  id_d  = NCDF_VARDEF(cdfid,'MODEL_INPUT_D',dim_i,/double)
  id_z0 = NCDF_VARDEF(cdfid,'MODEL_INPUT_Z0',dim_ii,/float) ; zero-plane displacement
  id_zm = NCDF_VARDEF(cdfid,'MODEL_INPUT_ZM',dim_ii,/float) ; effective measurement height
  id_L  = NCDF_VARDEF(cdfid,'MODEL_INPUT_OBUKHOV_LENGTH',dim_ii,/float) ; obukhov length
  id_sv = NCDF_VARDEF(cdfid,'MODEL_INPUT_SIGMA_V',dim_ii,/float) ; obukhov length
  id_u  = NCDF_VARDEF(cdfid,'MODEL_INPUT_WIND_VELOCITY',dim_ii,/float) ; obukhov length
  id_wd = NCDF_VARDEF(cdfid,'MODEL_INPUT_WIND_DIRECTION',dim_ii,/float) ; wind direction
   
  id_pm = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_M',dim_pp,/float) ; Exponent of the wind velocity power law
  id_pn = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_N',dim_pp,/float) ; Exponent of the eddy diffusivity power law
  id_uc = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_CONST_U',dim_pp,/float) ; Constant in power-law profile of the wind velocity
  id_us = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_USTAR',dim_pp,/float) ; Friction velocity
  id_ka = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_KAPPA',dim_pp,/float) ;Constant in power-law profile of the eddy diffusivity
  id_xi = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_XSI',dim_pp,/float) ;Flux length scale
  id_xm = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_XMAX',dim_pp,/float) ;Distance of maximum flux per unit point source 
  id_tp = NCDF_VARDEF(cdfid,'MODEL_PARAMETER_TOTALPHI',dim_pp,/float) ;Total of phi within whole grid
  
   ;define variable attributes --------------------------------
   
  NCDF_ATTPUT, cdfid, id_x,  'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_y,  'UNIT', 'metre', /char  
  NCDF_ATTPUT, cdfid, id_phi,'UNIT', 'metre^-2', /char 
  NCDF_ATTPUT, cdfid, id_phi,'DESCRIPTION', 'Cumulative flux footprint or vertical flux per unit point source'
  NCDF_ATTPUT, cdfid, id_d,  'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_d,  'DESCRIPTION', 'Spatial resolution of flux footprint raster data'
  NCDF_ATTPUT, cdfid, id_z0, 'UNIT', 'metre', /char
  NCDF_ATTPUT, cdfid, id_z0, 'DESCRIPTION', 'Roughness length (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_zm, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_zm, 'DESCRIPTION', 'Effective measurement height (z-d) (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_L , 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_L,  'DESCRIPTION', 'Obukhov length at zm (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Standard deviation of lateral wind component at zm (for each run in cumulative footprint)', /char
  NCDF_ATTPUT, cdfid, id_u , 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Longitudinal wind velocity component at zm (for each run in cumulative footprint)', /char
  NCDF_ATTPUT, cdfid, id_wd, 'UNIT', 'degree', /char
  NCDF_ATTPUT, cdfid, id_sv, 'DESCRIPTION', 'Wind direction relative to geographic North (clockwise) at zm (for each run in cumulative footprint)', /char
  NCDF_ATTPUT, cdfid, id_pm, 'UNIT', '-', /char  
  NCDF_ATTPUT, cdfid, id_pm, 'DESCRIPTION', 'Exponent of the wind velocity power law (for each run in cumulative footprint)', /char  
  NCDF_ATTPUT, cdfid, id_pn, 'UNIT', '-', /char 
  NCDF_ATTPUT, cdfid, id_pn, 'DESCRIPTION', 'Exponent of the eddy diffusivity power law (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_uc, 'UNIT', 'metre^(1−m) second^−1', /char 
  NCDF_ATTPUT, cdfid, id_uc, 'DESCRIPTION', 'Constant in power-law profile of the wind velocity (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_us, 'UNIT', 'metre second^-1', /char 
  NCDF_ATTPUT, cdfid, id_us, 'DESCRIPTION', 'Friction velocity', /char 
  NCDF_ATTPUT, cdfid, id_ka, 'UNIT', 'metre^(2−n) second^−1', /char 
  NCDF_ATTPUT, cdfid, id_ka, 'DESCRIPTION', 'Constant in power-law profile of the eddy diffusivity (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_xi, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_xi, 'DESCRIPTION', 'Flux length scale', /char 
  NCDF_ATTPUT, cdfid, id_xm, 'UNIT', 'metre', /char 
  NCDF_ATTPUT, cdfid, id_xm, 'DESCRIPTION', 'Upwind distance of maximum flux per unit point source (for each run in cumulative footprint)', /char 
  NCDF_ATTPUT, cdfid, id_tp, 'UNIT', 'domain^-2', /char 
  NCDF_ATTPUT, cdfid, id_tp, 'DESCRIPTION', 'Total phi within the boundaries of the grid (for each run in cumulative footprint)', /char 
  
  ;put the file in data mode ---------------------------------

  NCDF_CONTROL, cdfid, /endef

  ;put data in variable --------------------------------------
  
  NCDF_VARPUT, cdfid, id_x, easting
  NCDF_VARPUT, cdfid, id_y, northing
  NCDF_VARPUT, cdfid, id_phi, cummulative_phi
  NCDF_VARPUT, cdfid, id_d, d_input
  NCDF_VARPUT, cdfid, id_z0, z_0_input
  NCDF_VARPUT, cdfid, id_zm, z_m_input
  NCDF_VARPUT, cdfid, id_L, L_input
  NCDF_VARPUT, cdfid, id_sv, sigv_input
  NCDF_VARPUT, cdfid, id_u, u_input
  NCDF_VARPUT, cdfid, id_wd, wd_input
  
  NCDF_VARPUT, cdfid, id_pm, param_m
  NCDF_VARPUT, cdfid, id_pn, param_n
  NCDF_VARPUT, cdfid, id_uc, param_u
  NCDF_VARPUT, cdfid, id_us, ustar
  NCDF_VARPUT, cdfid, id_ka, param_kappa
  NCDF_VARPUT, cdfid, id_xi, xsi
  NCDF_VARPUT, cdfid, id_xm, x_max
  NCDF_VARPUT, cdfid, id_tp, total_phi
  
  ;close ncdf file -------------------------------------------

  NCDF_CLOSE, cdfid
  
 
end


;+
; name:
;   ubc_fil_ncdf2struct.pro
;
; purpose:
;   reads variables, parameter attributes and global attributes from
;   ncdf-files into a single idl-structure.
;
; category:
;   file handling
;
; calling sequence:
;   ubc_fil_ncdf2struct, ncdf_file
;
; input:
;   ncdf_file : string. full filepath of ncdf file.
;   data : structure according to ncdf-file.
;
; restrictions:
;   dimension-names are not retreived.
;
; examples
;   data=ubc_fil_ncdf2struct('c:\temp\test.ncdf')
;   help, data, /stru
;
; modification history:
;   andi christen, tu berlin, 23-aug-05, andreas.christen@tu-berlin.de
;-

function ncdf2struct, ncdf_file

  ;open ncdf-file -----------------------------------------

  cdfid = ncdf_open(ncdf_file)

  ;inquire ncdf-file --------------------------------------

  inq = ncdf_inquire(cdfid)

  ;resolve parameter names ---------------------------------

  for varid=0, inq.nvars-1 do begin
    varinq = ncdf_varinq(cdfid, varid)
    if varid eq 0 then begin
      var_name=varinq.name
    endif else var_name=[var_name,varinq.name]
  endfor

  ;resolve parameter values and parameter attributes --------

  for varid=0, inq.nvars-1 do begin
    varinq = ncdf_varinq(cdfid, varid)
    ncdf_varget, cdfid, varid, value
    if strlowcase(varinq.datatype) eq 'char' then value=string(value)
    if varid eq 0 then variables = create_struct(var_name[varid],value) else $
      variables = create_struct(variables,var_name[varid],value)
    for va=0, varinq.natts-1 do begin
      vattname = ncdf_attname(cdfid, varid, va)
      var_attb_inq = ncdf_attinq(cdfid, varid, vattname)
      ncdf_attget, cdfid, varid, vattname, value
      if strlowcase(var_attb_inq.datatype) eq 'char' then value=string(value)
      variables = create_struct(variables,var_name[varid]+'_'+vattname,value)
    endfor
  endfor

  ;resolve global attributes ------------------------------

  attributes={name:'',datatype:'',length:0l}
  if inq.ngatts gt 0 then begin
    attributes=replicate(attributes,inq.ngatts)
    for a=0, inq.ngatts-1 do begin
      attn = ncdf_attname(cdfid, /global , a)
      attributes[a].name=attn
      att_stru=ncdf_attinq(cdfid,/global,attn)
      attributes[a].datatype=strupcase(att_stru.datatype)
      attributes[a].length=att_stru.length
    endfor

    for a=0,inq.ngatts-1 do begin
      ncdf_attget, cdfid, /global, attributes[a].name, value
      if strlowcase(attributes[a].datatype) eq 'char' then value=string(value)
      if a eq 0 then global_attributes = create_struct(attributes[a].name,value) else $
        global_attributes = create_struct(global_attributes,attributes[a].name,value)
    endfor
  endif

  ;close ncdf file -----------------------------------------

  ncdf_close, cdfid

  ;merge variables / attributes into single structure ------

  if inq.ngatts gt 0 then data=create_struct(global_attributes,variables) else $
    data=variables

  return, data

end
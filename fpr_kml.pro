;+
; name: fpr_kml
;
; purpose:
;   translates a netCDF footprint file into a .kml file that can be opened in Google Earth
;   for visualization purposes.
;
; category:
;   Footprint modelling
;
; calling sequence:
;   ubc_fpr_kml, ncfile, lon, lat
;
; required inputs:
;   ncfile     : path to a single gridded netCDF footprint file created with fpr_write_ncdf.pro 
;   (can be a cumulative file)
;   lon = longitude of EC-system / tower location in WGS-84 Lat/Lon (West = negative, East = positive)
;   lat = latitude of EC-system / tower location in WGS-84 Lat/Lon (South = negative, North = positive)
;   
; output:
;   .kml file and .png file as layer to be included. .kml file and .png file need to be in the same
;   directory for proper visualization.
;
; example:
;  ;individual calculation
;  fpr_kml, '/Users/name/Desktop/test.nc', -123.0784, 49.2261, /heat_map
;  
; dependencies
;   needs coyote IDL library for "cgUTMZone" and "cgMap"
;   needs ncdf2struct.pro from https://github.com/achristen/netCDF2struct-IDL
; 
; revision history:
;   2016-12-01 ac
;-


pro fpr_kml, ncfile, lon, lat, heat_map=heat_map
 
  if keyword_set(heat_map) then viz_type = '_HM' else viz_type=''
  
  kmlfile = FILE_DIRNAME(ncfile)+path_sep()+file_basename(ncfile,'.nc')+viz_type+'.kml'
  pngfile = FILE_DIRNAME(ncfile)+path_sep()+file_basename(ncfile,'.nc')+viz_type+'.png'
  
  if n_params() lt 1 then ncfiles=dialog_pickfile(default_extension='nc',/read,title='Select Footprint NCDF...')
  footprint = ncdf2struct(ncfile[0])

  so = reverse(sort(footprint.phi))
  summed = footprint.phi & summed[*] = 0D
  nx = n_elements(footprint.easting)
  ny = n_elements(footprint.northing)
  
  sum = 0D
  
  for s=0L, n_elements(so)-1 do begin
    sum = sum + footprint.phi[so[s]] 
    summed[so[s]] = sum
  endfor
  
  Description = footprint.description+'. Individual time steps included: '+strcompress(string(n_elements(footprint.juliantime),format='(i10)'),/remove_all)+'. Model resolution: d = '+strcompress(string(footprint.MODEL_INPUT_D,format='(f10.1)'),/remove_all)+'m'
  
  min_jul=min(footprint.juliantime)
  max_jul=max(footprint.juliantime)
  caldat, min_jul, min_day, min_mon, min_year, min_hour, min_min
  caldat, max_jul, max_day, max_mon, max_year, max_hour, max_min
  min_sdat = strupcase(string(min_year,format='(i4.4)')+'-'+string(min_mon,format='(i2.2)')+'-'+string(min_day,format='(i2.2)')+' '+string(min_hour,format='(i2.2)')+':'+string(min_min,format='(i2.2)'))
  max_sdat = strupcase(string(max_year,format='(i4.4)')+'-'+string(max_mon,format='(i2.2)')+'-'+string(max_day,format='(i2.2)')+' '+string(max_hour,format='(i2.2)')+':'+string(max_min,format='(i2.2)'))

  if min_jul eq max_jul then time_str = min_sdat else $
  time_str = strupcase(min_sdat + ' to ' + max_sdat )
     
  contours = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.90,0.95] 
  color_set = ['medium violet red','indigo','dark blue','royal blue','teal','green','lime green','yellow green','pale goldenrod','white','white','white','white']
  
  ;get border pixels to determine if any countour is outside bounds of grid
  border_index = bytarr(nx,ny)
  border_index[*] = 0B
  border_index[0,*] = 1B
  border_index[*,0] = 1B
  border_index[nx-1,*] = 1B
  border_index[*,ny-1] = 1B
  border_loc = where(border_index eq 1)
  
  highest_percentage_at_border = min(summed[border_loc])  
  valid_contour_index = where(contours lt highest_percentage_at_border, cnt)
  contours = contours[valid_contour_index]
  color_set = color_set[valid_contour_index]
  color_set[n_elements(color_set)-1] = 'white'
  
  if keyword_set(heat_map) then begin
    
  lowest_alog = -11
  highest_alog = -7
    
  ; heat map of source area
  
  transformed = alog(footprint.phi)
  below = where(transformed lt lowest_alog, cnt)
  if cnt gt 0 then transformed[below] = lowest_alog
  above = where(transformed gt highest_alog, cnt)
  if cnt gt 0 then transformed[above] = highest_alog
  transformed = byte((254.0/(highest_alog-lowest_alog))*(transformed-lowest_alog))
  
  transparency = transformed
  a = image(transformed, footprint.easting, footprint.northing, margin=[0.0,0.0,0.0,0.0], rgb_table=0, dimension=[1000,1000], /buffer) 
  op = plot([footprint.easting[0],footprint.easting[0],footprint.easting[nx-1],footprint.easting[nx-1],footprint.easting[0]]*0.998,$
    [footprint.northing[0],footprint.easting[ny-1],footprint.easting[ny-1],footprint.easting[0],footprint.easting[0]]*0.998,AXIS_STYLE=0,color=255,thick=4,linestyle=2,/overplot)

  a.save, pngfile, RESOLUTION = 180, BIT_DEPTH=1
  a.close
  
  back_in = read_png(pngfile, R, G, B)
  transpmask = r
  R = bytarr(n_elements(r)) + 255B
  G = bytarr(n_elements(r)) + 0B
  B = bytarr(n_elements(r)) + 0B
  
  write_png, pngfile, back_in, R, G, B, transparent=transpmask
  
  endif else begin

  ; traditional contour plot of source area

  a = contour(summed, footprint.easting, footprint.northing, margin=[0.0,0.0,0.0,0.0], dimension=[1000,1000], $
    /fill, xticklen=-0.01, yticklen=-0.01, ASPECT_RATIO=1.0, BACKGROUND_TRANSPARENCY=100, $
    xtitle='Easting (m)', ytitle='Northing (m)', c_value=contours, C_LABEL_SHOW=0, BACKGROUND_COLOR='white', $
    C_COLOR=color_set, /buffer)
  a = contour(summed*100, footprint.easting, footprint.northing, LABEL_FORMAT='(i2)', c_value=contours*100, color='black',/overplot)
  
  op = plot([footprint.easting[0],footprint.easting[0],footprint.easting[nx-1],footprint.easting[nx-1],footprint.easting[0]]*0.998,$
     [footprint.northing[0],footprint.easting[ny-1],footprint.easting[ny-1],footprint.easting[0],footprint.easting[0]]*0.998,AXIS_STYLE=0,color='black',thick=4,linestyle=2,/overplot)

  a.save, pngfile, RESOLUTION = 180, BIT_DEPTH=1
  a.close
  
  back_in = read_png(pngfile, R, G, B)
  white = where(r eq 248 and g eq 248 and b eq 248)
  transpmask = bytarr(n_elements(r))+255B
  transpmask[white] = 0B
  write_png, pngfile, back_in, R, G, B, transparent=transpmask
  
  endelse
  
  mapCoord = cgMap('UTM', Zone=cgUTMZone(lon,lat), Ellipsoid='WGS84')
  centre = mapCoord.Forward(lon,lat)
  
  d_half = footprint.model_input_d / 2
  sw_corner = [centre[0]+footprint.easting[0]-d_half, centre[1]+footprint.northing[0]-d_half]
  se_corner = [centre[0]+footprint.easting[n_elements(footprint.easting)-1]+d_half, centre[1]+footprint.northing[0]-d_half]
  ne_corner = [centre[0]+footprint.easting[n_elements(footprint.easting)-1]+d_half, centre[1]+footprint.northing[n_elements(footprint.northing)-1]+d_half]
  nw_corner = [centre[0]+footprint.easting[0]-d_half, centre[1]+footprint.northing[n_elements(footprint.northing)-1]+d_half]
  
  sw_corner_wgs84 = mapCoord.Inverse(sw_corner[0], sw_corner[1])
  se_corner_wgs84 = mapCoord.Inverse(se_corner[0], se_corner[1])
  ne_corner_wgs84 = mapCoord.Inverse(ne_corner[0], ne_corner[1])
  nw_corner_wgs84 = mapCoord.Inverse(nw_corner[0], nw_corner[1])
  
  print, sw_corner_wgs84
  print, se_corner_wgs84
  print, ne_corner_wgs84
  print, nw_corner_wgs84
  
  openw, lun, kmlfile, /get_lun
  printf, lun, '<?xml version="1.0" encoding="UTF-8"?>'
  printf, lun, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">'
  printf, lun, '<GroundOverlay>'
  printf, lun, '<name>'+file_basename(ncfile,'.nc')+'</name>'
  printf, lun, '<description>'+Description+'</description>'
  if not keyword_set(heatmap) then printf, lun, '<color>80ffffff</color>' ; make transparent partially
  printf, lun, '<Icon>'
  printf, lun, '<name>Footprint</name>'
  printf, lun, '  <href>'+file_basename(pngfile)+'</href>'
  printf, lun, '  <viewBoundScale>0.75</viewBoundScale>'
  printf, lun, '</Icon>'
  printf, lun, '<gx:LatLonQuad>'
  printf, lun, '  <coordinates>'
  printf, lun, strcompress(string(sw_corner_wgs84[0])+','+string(sw_corner_wgs84[1])+',0',/remove_all)+' '+$
               strcompress(string(se_corner_wgs84[0])+','+string(se_corner_wgs84[1])+',0',/remove_all)+' '+$
               strcompress(string(ne_corner_wgs84[0])+','+string(ne_corner_wgs84[1])+',0',/remove_all)+' '+$
               strcompress(string(nw_corner_wgs84[0])+','+string(nw_corner_wgs84[1])+',0',/remove_all)
  printf, lun, '  </coordinates>'
  printf, lun, '</gx:LatLonQuad>'
  printf, lun, '</GroundOverlay>'
  printf, lun, '</kml>'
  close, lun
  free_lun, lun
  
end
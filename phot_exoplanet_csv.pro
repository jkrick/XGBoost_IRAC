pro phot_exoplanet_csv, planetname, apradius,chname, rawfile = rawfile
;
;+
; PROCEDURE:
;    phot_exoplanet_csv
;
; PURPOSE: 
;    Do photometry on any IRAC warm mission staring mode exoplanet data
;
; USAGE:
;      phot_exoplanet_csv, 'XO3', 2.25, '2',/rawfile
;
; INPUT:
;    planetname: strong containing the name of the planet
;      apradius: aperture radius at which to do the photometry [ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
;     chname: string containing the number of the IRAC channel: '1' or '2'
;
; OUTPUTS:
;    csv file with requested photometry
;
; OPTIONAL FLAGS:
;    rawfile: if set, use warm raw data to determine pixel values
;    (instead of bcd)
;
; COMMENTS:
;    A lot of care is taken to make the distinction between sub array
;    and full array data.  Sub aray has naxis = 3 set in the header
;    because it is a data cube of 64 images which are 32x32.  Full
;    array has naxis = 2 set in the header because it is just a
;    regular image.
  
 
basedir = '/Users/jkrick/external/irac_warm/' 
aorname= ['r46467072', 'r46471424', 'r46467840', 'r46471168', 'r46470144', 'r46470912', 'r46467584', 'r46470656', 'r46469376','r46470400', 'r46466816', 'r46469632', 'r46468864', 'r46469120', 'r46469888', 'r46468608', 'r46467328', 'r46468352','r46471680', 'r46468096']

;setup

;convert aperture radius in pixels into what
;the aperture photometry routine uses
case apradius of
   ;[ 2.0, 2.25, 2.5, 2.75, 3.0]
   2.0: apval = 0
   2.25: apval = 1
   2.5:  apval = 2
   2.75: apval =3
   3.0: apval = 4
   
   Else: apval = 2              ; if no decision, then choose an apradius = 2.5 pixels
endcase

;;choose pmap file based on channel for use in correcting pixel phase effect
if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_rmulti_s3_7_0p1s_x4_150723.sav'
if chname eq '1' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_rmulti_s3_7_sd_0p4s_sdark_150722.sav'


;---------------

dirname = strcompress(basedir + planetname +'/')

startaor =0
stopaor =   n_elements(aorname) - 1
for a =startaor,  stopaor do begin  ; for each AOR
   print, '----------------------------'
   print, 'working on ',a, ' ', aorname(a)
   dir = dirname+ string(aorname(a) ) 
   CD, dir                      ; change directories to the correct AOR directory

   ;;find all the fits files and make a list of each type of file
   command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
   spawn, command
   command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
   spawn, command2
   command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
   spawn, command3

   ;;read in the file lists
   readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
   readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
   readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent

   print,'n_elements(fitsname)', n_elements(fitsname)   

   startfits = 0L
   for i =startfits,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
      
      ;;start by just reading in the header and setting up some arrays
      ;;for bookkeeping.
      header = headfits(fitsname(i)) ;
      sclk_obs= sxpar(header, 'SCLK_OBS')
      frametime = sxpar(header, 'FRAMTIME')
      bmjd_obs = sxpar(header, 'BMJD_OBS')
      ch = sxpar(header, 'CHNLNUM')
      ronoise = sxpar(header, 'RONOISE') ; these are zeros
      gain = sxpar(header, 'GAIN')
      fluxconv = sxpar(header, 'FLUXCONV')
      exptime = sxpar(header, 'EXPTIME')
      aintbeg = sxpar(header, 'AINTBEG')
      atimeend = sxpar(header, 'ATIMEEND')
      afpat2b = sxpar(header, 'AFPAT2B')
      afpat2e = sxpar(header, 'AFPAT2E')
      afpat2 = (afpat2b + afpat2e)/2.      
      naxis = sxpar(header, 'NAXIS')
      framedly = sxpar(header, 'FRAMEDLY')
      ra_ref = sxpar(header, 'RA_REF')
      dec_ref = sxpar(header, 'DEC_REF')

      if ch eq '2' and frametime eq 2 then ronoise = 12.1
      if i eq startfits then sclk_0 = sclk_obs
      
      if i eq startfits and naxis eq 3 then begin
         ;;still setting up bookkeeping
         xarr = fltarr(63*(n_elements(fitsname)))
         yarr = xarr
         xerr = xarr
         yerr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
         corrfluxarr = xarr
         corrfluxerrarr = xarr
         timearr =  dblarr(63*(n_elements(fitsname)))
         bmjd = dblarr(63*(n_elements(fitsname)))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         temparr = xarr
         edgemean = xarr
         piarr = fltarr(7,7,63L*n_elements(fitsname))
      endif
      if i eq startfits and naxis ne 3 then begin
         xarr = fltarr(n_elements(fitsname))
         yarr = xarr
         xerr = xarr
         yerr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
         corrfluxarr = xarr
         corrfluxerrarr = xarr
         timearr =dblarr(n_elements(fitsname))
         bmjd = dblarr(n_elements(fitsname))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         temparr = xarr
         piarr =fltarr(7,7,n_elements(fitsname))
      endif
      fdarr = fltarr(n_elements(fitsname))
      fdarr[i] = framedly

      if naxis eq 3 then begin
         deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
         nt = dindgen(64)
         sclkarr = sclk_obs  + (deltatime*nt)
         bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D  
      endif else begin          ; full array, so don't need to expand out the times
         sclkarr = sclk_obs
         bmjdarr = bmjd_obs
      endelse

      ;read in the fits files
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc      

      ;apply mask file if necessary to mask out neighboring stars
      if planetname eq 'hd189733' then  im[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
      if planetname eq 'HAT-P-22' and chname eq '2' then  im[19:25, 9:15, *] = !Values.F_NAN  
      if planetname eq 'HAT-P-22' and chname eq '1' then  im[5:11, 11:17, *] = !Values.F_NAN  
      if planetname eq 'HD93385' then im[17:21, 21:27, *] = !Values.F_NAN                      
      
      if planetname eq 'WASP-14b' then  begin
         ra_off = 218.27655     ; 218.27718
         dec_off = 21.897652    ; 21.89808
         adxy, h, ra_off, dec_off, x_off, y_off
         x_off = round(x_off)
         y_off = round(y_off)
         im[(x_off-3):(x_off + 3), (y_off - 3):(y_off+3), *] = !Values.F_NAN ;mask region with nan set for bad regions
      endif
      
      ;run the centroiding and photometry
      get_centroids,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm,  /WARM
      nanfound = where(FINITE(np) lt 1, nancount)
      if nancount gt 0 then print, 'NAN: ', fitsname(i), nanfound
      
      x_center = temporary(x3)
      y_center = temporary(y3)
      xcenerr = temporary(x3s)
      ycenerr = temporary(y3s)
      
     ;choose the requested pixel aperture
      abcdflux = f[*,apval]      
      fs = fs[*,apval]
     ; 3-7 pixel background
      back = b[*,0]
      backerr = bs[*,0]
      npcentroids = np         ;noise pixel

      
      if keyword_set(rawfile) then begin
         ;;read in the raw data file and get the DN level of the peakpixel      
         fits_read, rawname(i), rawdata, rawheader
        
         barrel = fxpar(rawheader, 'A0741D00')
         fowlnum = fxpar(rawheader, 'A0614D00')
         pedsig = fxpar(rawheader, 'A0742D00')
         ichan = fxpar (rawheader, 'CHNLNUM')
         tb = fxpar(rawheader, 'A0655D00')
         te = fxpar(rawheader, 'A0664D00')
         mean_t = (tb + te)/2.
         
         ;;use Bill's code to convert to DN
         dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
         rawdata = reform(rawdata, 32, 32, 64)
         peakpixDN = abcdflux   ; set up the array to be the same size as abcdflux
         if naxis eq 3 then begin
            for pp = 0, 63 do begin
               peakpixDN[pp] = max(rawdata[13:16,13:16,pp])
            endfor
         endif

         ;;convert temperature DN into degrees kelvin
         t_cernox = ch EQ '1' ? fpa1t2(mean_t,4) : fpa2t2(mean_t,4)

      endif

      ;;track the values of the 7x7 pixel box around the centroid
      if naxis eq 3 then begin
         pi = im[12:18, 12:18,*] ; now a 7x7x64 element array
         edgepix = [im[12,12,*], im[12,13,*],im[12,14,*],im[12,15,*],im[12,16,*],im[12,17,*],im[12,18,*],im[18,12,*], im[18,13,*],im[18,14,*],im[18,15,*],im[18,16,*],im[18,17,*],im[18,18,*],im[13,12,*],im[14,12,*], im[15,12,*], im[16,12,*], im[17,12,*],im[13,18,*],im[14,18,*], im[15,18,*], im[16,18,*], im[17,18,*]]
         edge = reform(edgepix)
        edge_mean =  mean(edge,dimension = 1, /nan)
      endif
      
      if naxis eq 2 then begin
         pi = im[(fix(x_center) - 3):(fix(x_center+3)), (fix(y_center) - 3):(fix(y_center+3))] ; now a 7x7  array

      endif
      
   
;---------------------------------
      ;;use pmap data to find nearest neighbors in pmap dataset and find a
      ;;correction based on those neighbors.  this is referred to as the
      ;;pixel phase correction.
      if naxis eq 3 then begin
         corrflux = pmap_correct(x_center,y_center,abcdflux,ch,xfwhm,yfwhm,NP = npcentroids,$
                                 FUNC=fs,CORR_UNC=corrfluxerr, DATAFILE=pmapfile,NNEAREST=nn, $
                                 R_USE = apradius, USE_PMAP = IMAIN) ;,/VERBOSE
      endif else begin
         corrflux = pmap_correct(x_center,y_center,abcdflux,ch,xfwhm,yfwhm,NP = npcentroids,$
                                 FUNC=fs,CORR_UNC=corrfluxerr, DATAFILE=pmapfile,NNEAREST=nn, $
                                 R_USE = apradius, USE_PMAP = IMAIN,/full) ;,/VERBOSE
      endelse
         


;---------------------------------
      ;;if subarray ignore the first image which is known to have some problems
      if naxis eq 3 then begin 
         xarr[i*63] = x_center[1:*]
         yarr[i*63] = y_center[1:*]
         xerr[i*63] = xcenerr[1:*]
         yerr[i*63] = ycenerr[1:*]
         fluxarr[i*63] = abcdflux[1:*]
         fluxerrarr[i*63] = fs[1:*]
         corrfluxarr[i*63] = corrflux[1:*]
         corrfluxerrarr[i*63] = corrfluxerr[1:*]
         timearr[i*63] = sclkarr[1:*]        
         bmjd[i*63] = bmjdarr[1:*]
         backarr[i*63] = back[1:*]
         backerrarr[i*63] = backerr[1:*]
         nparr[i*63] = np[1:*]
         npcentroidsarr[i*63] = npcentroids[1:*]
         xfwhmarr[i*63] = xfwhm[1:*]
         yfwhmarr[i*63] = yfwhm[1:*]
         temparr[i*63] =fltarr(63) + t_cernox        ;;make an array with 63 elements of the temperature
         edgemean[i*63] = edge_mean[1:*]
         
         if keyword_set(rawfile) then peakpixDNarr[i*63] = peakpixDN[1:*]

         ;;need to loop through each image individually
         for pj = 0, 62 do begin
            piarr(*,*,i*63 + pj) = pi[*,*,pj]
         endfor
         
 
      endif 
      if naxis eq 2 then begin  ; and i eq 0 then begin
         xarr[i] = x_center
         yarr[i]  =  y_center
         yerr[i] = ycenerr
         xerr[i] = xcenerr
         fluxarr[i]  =  abcdflux
         fluxerrarr[i]  =  fs
         corrfluxarr[i]  = corrflux
         corrfluxerrarr[i]  =  corrfluxerr
         timearr[i]  = sclkarr
         bmjd[i]  = bmjdarr
         backarr[i]  =  back
         backerrarr[i]  = backerr
         nparr[i]  = npcentroids
         npcentroidsarr[i] = npcentroids
         xfwhmarr[i] = xfwhm
         yfwhmarr[i] = yfwhm
         temparr[i] = t_cernox
         
         if keyword_set(rawfile) then peakpixDNarr[i] = peakpixDN
         piarr[*,*,i] = pi
      endif

      if a eq startaor and i eq startfits then begin
         time_0 = bmjdarr(0)
      endif

   endfor; for each fits file in the AOR


   
  
;--------------------------------
;write out csv file
;--------------------------------
   ;;make a 2D array to hold variables to output to csv file
   pix2_2 = reform(piarr[2, 2, *])
   pix2_3 = reform(piarr[2, 3, *])
   pix2_4 = reform(piarr[2, 4, *])
   pix3_2 = reform(piarr[3, 2, *])
   pix3_3 = reform(piarr[3, 3, *])
   pix3_4 = reform(piarr[3, 4, *])
   pix4_2 = reform(piarr[4, 2, *])
   pix4_3 = reform(piarr[4, 3, *])
   pix4_4 = reform(piarr[4, 4, *])
   
   
   all_arr = [[xarr],[xerr], [yarr],[yerr], [fluxarr], [fluxerrarr], [npcentroidsarr], [xfwhmarr], [yfwhmarr],  [bmjd],  [backarr], [backerrarr], [pix2_2] , [pix2_3], [pix2_4], [pix3_2] , [pix3_3], [pix3_4], [pix4_2] , [pix4_3], [pix4_4],[temparr], [edgemean]]
    
   trans_arr = transpose(all_arr)

   header = ['xpos','xerr','ypos','yerr','flux', 'fluxerr','np', 'xfwhm','yfwhm', 'bmjd',  'bg_flux', 'sigma_bg_flux','pix1','pix2', 'pix3','pix4','pix5','pix6','pix7','pix8','pix9','t_cernox','edge_mean']
   
   ;;use write_csv
   csvname = strcompress('/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs'+ '/'+ planetname + '/' + planetname +'_'+ aorname(a)+'.csv')
   write_csv, csvname, trans_arr, header = header
   
    
endfor                          ;for each AOR



 
end




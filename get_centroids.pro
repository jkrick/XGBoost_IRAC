pro get_centroids, im, h, unc, ra, dec, t, dt, hjd, xft, x3, y3, $
                       x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                       x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                       xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                       xfwhmarr, yfwhmarr, bb, WARM=warm, SILENT=silent
                       
;+
; PROCEDURE:
;    GET_CENTROIDS
;
; PURPOSE: 
;    Return x and y position of brightest source for an IRAC
;    image.  Will return all instances of source if a subarray frame is
;    supplied.  Should also work on pbcd mosaics and other files.  This version
;    does not read in files, but needs the data passed to it.
;
; USAGE:
;      fits_read, fitsname(i), im, h
;      fits_read, buncname(i), unc, hunc
;      get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
;                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
;                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
;                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
;                                   xfwhm, yfwhm, /WARM
; INPUT:
;    im: array containing  data
;      h: header of the data fits file
;     unc: array containing the uncertainty image data
;
; OUTPUTS:
;    t: double array containing SCLK times of subarray images
;    dt: double array containing time between successive images
;    x0: float array containing x coordinate of centroid
;    y0: float array containing y coordinate of centroid
;     f: 2D float array containing flux of source in specified aperture
 ; f[i,1] = fluxes in all subarray frames with [3,12,20] aperture photometry
;    xs: float array containing uncertainty in x coordinate of centroid
;    ys: float array containing uncertainty in y coordinate of centroid
;    fs: float array containing flux uncertainty of source
;     b: float array containing background
;    bs: float array containing uncertainty in background
;     c: array containing number of pixels in source aperture
;    cb: array containing number of pixels in background aperture
;
; OPTIONAL FLAGS:
;    WARM: if set, use warm IRAC factors for noise
;
; METHOD:
;    The brightest pixel in each slice of the data cube is found, a median
;    background is removed and the first moments are calculated using a 
;    2*boxwidth+1 box centered on the brightest pixel.  aper.pro is also used if 
;    desired to perform circular aperture photometry
;    The apertures used are [2, 3, 4, 5, 6, 8, 10, 12] 
;    with a [12, 20] background annulus, 2 with a [2, 6] background, 3 with 
;    a [3, 7] background, 5 with a [5, 10] background.  Box centroiding using 3, 5, 7
;    pixel boxes for centroid measures as well as gcntrd and cntrd from IDLphot.

;
; HISTORY:
;    27 July 2011 Added fwhm
;    11 May 2011 SJC Hard code full array plate scale as those are the ones we use for both
;                full and subarray photometry
;    15 Apr 2011 SJC Increased saturation levels by 5%
;    12 Dec 2009 SJC Modified to not read files, but handle input images
;    24 Nov 2009 SJC Added code so that pBCD mosaics are handled, but still need 
;                to appropriately use coverage map
;    03 August 2009 SJC Modified for warm mission
;    17 Jan 2008 SJC Corrected typos in code. 
;    13 Nov 2006 SJC Added flux an03 aud flux uncertainty measurement, added 
;                keywords for background annulus and box border
;    14 September 2004 SJC Added uncertainty calculations
;    08 September 2004 SJC Cleanup of initial code
;    10 March 2005 SJC Corrected calculation of sclk for raw data
;    June 2004 SJC Initial code 
;-
	name = 'GET_CENTROIDS_FOR_CALSTAR:'

; GOO Syntax is out of date
;	str = 'Syntax - get_centroids_for_calstar, file, t, dt, hjd, x0, y0, [f, xs, ys, $';
;	str1 = '                   fs,b,bs,c,cb,cube], $'
;	str2 = '                   SIGFILE=sigfile, $'
;	str3 = '                   SUPSIGMA=supsigma, /WARM, $'
;	str4 = '                   /APER, /SILENT, /ONLY_POS, RA=ra, DEC=dec, $'
;	str5 = '                   NOISE_PIXELS=np, COVERAGE=cfile, /BIG'
;	str6 = 'Syntax is out of date'
;	if (N_params() lt 6 or N_params() gt 15) then begin
;		print, str
;		print, str1
;		print, str2
;		print, str3
;		print, str4
;		print, str5		
;		return
;	endif

	if (not keyword_set(SILENT)) then silent = 0 else silent = 1

; Bad pixel settings needed for aper.pro
	badpix = [-9., 9.] * 1.D8
	
; Edge limit
	edge = 5.

; First set of apertures	
  aps1 = [ 2.0, 2.25, 2.5, 2.75, 3.0];, 2.50, 3.0];, 2.5];, 2.75]
;  aps1 = [ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
  naps1 = n_elements(aps1) 
; First background annulus
	back1 =  [3., 7.]
; Second set of apertures
	aps2 = [2.25]
; Second background annulus
	back2 = [3.,7.]
; Third set of apertures
	aps3 = [3.]
; Third background annulus
	back3 = [3., 7.]
; Fourth set of apertures
	aps4 = [5.]
; Fourth background annulus
	back4 = [5., 10.]
        back5 = [7., 14.]
        back6 = [3., 15.]
; Number of apertures
	napers = n_elements(aps1) ;+ 2
; Number of background annuli
	nbacks = 1; 3;4

; full width at half maximum of sources imaged by IRAC in arcseconds
	fwhm = [1.66, 1.72, 1.88, 1.98]
	
; Plate scales for full array BCDs
	pxscal1 = [-1.22334117768332D, -1.21641835430637D, -1.22673962032422D, -1.2244968325831D]
	pxscal1 = abs(pxscal1)
	pxscal2 = [1.22328355209902D, 1.21585676679388D, 1.22298117494211D, 1.21904995758086D]


; Array of readnoise in electrons, per channel for 0.02, 0.1, 0.4 frametimes
; 2s sub, 0.4, 0.6, 1.2, 2, 6, 12, 30, 100 second full
	if (keyword_set(WARM)) then begin
		readnoise = [[23.9, 23.8, 12.7, 11.9], [16.9, 16.8, 9.0, 8.4], $
	             [11.8, 12.1, 9.1, 7.1], [9.4, 9.4, 8.8, 6.7], $
	             [23.9, 23.8, 12.7, 11.9], [23.9, 23.8, 12.7, 11.9], $
	             [23.9, 23.8, 12.7, 11.9], [11.8, 12.1, 9.1, 7.1], $
	             [9.4, 9.4, 8.8, 6.7], [9.4, 9.4, 8.8, 6.7], $
	             [6.0, 8.4, 4.5, 4.2], [6.0, 8.4, 4.5, 4.2]] * 1.D 
; Bump up readnoise by 5% for warm mission for FPA 1
	    readnoise[0, *] = readnoise[0, *] * 1.05	    
	endif else $
		readnoise = [[23.9, 23.8, 12.7, 11.9], [16.9, 16.8, 9.0, 8.4], $
	             [11.8, 12.1, 9.1, 7.1], [9.4, 9.4, 8.8, 6.7], $
	             [23.9, 23.8, 12.7, 11.9], [23.9, 23.8, 12.7, 11.9], $
	             [23.9, 23.8, 12.7, 11.9], [11.8, 12.1, 9.1, 7.1], $
	             [9.4, 9.4, 8.8, 6.7], [9.4, 9.4, 8.8, 6.7], $
	             [8.5, 8.4, 4.5, 4.2], [8.5, 8.4, 4.5, 4.2]] * 1.D

; Array of saturation limits per channel for all frametimes in mJy
; mm is limit for inappropriate frametime for channel or mission status
	mm = 0.001
; 0.02, 0.1, 0.4, 2.0 sub, 0.4, 0.6, 1.2, 2, 6, 12, 30, 100 frametimes
	if (keyword_set(WARM)) then $
		sat_levels = [[20000., 28000., mm, mm], [3000., 3700., mm, mm], $
	             [700., 900., mm, mm], [150., 180., mm, mm], $
	             [900., 1000., mm, mm], [500., 600., mm, mm], $
	             [250., 300., mm, mm], [160., 200., mm, mm], $
	             [50., 60., mm, mm], [25., 30., mm, mm], $
	             [10., 12., mm, mm], [2.8, 3.5, mm, mm]] $
	else $
		sat_levels = [[20000., 17000., 63000., 45000.], [4000., 3300., 13000., 9000.], $
	             [1000., 820., 3100., 2300.], [mm, mm, mm, mm], $
	             [950., 980., 6950., 3700.], [630., 650., 4600., 2500.], $
	             [23.9, 23.8, 12.7, 11.9], [190., 200., 1400., 740.], $
	             [mm, mm, mm, mm], [32., 33., 230., 120.], $
	             [13., 13., 92., 48.], [3.8, 3.9, 27., 28.]]	
	             
; Increase saturation levels by 5% to account for them being pessimistic
	sat_levels = sat_levels * 1.05
; Gain (conversion from DN to e-) per channel
	if (keyword_set(WARM)) then gain = [3.7, 3.7, 3.8, 3.8] * 1.D$
	else gain = [3.3, 3.7, 3.8, 3.8] * 1.D
	
; Sigma level for outlier rejection
	sigma = 3.0


; Get channel number of image
	nch = sxpar(h, 'CHNLNUM')

; Set image dimensions
        nx = sxpar(h, 'NAXIS1')
        ny = sxpar(h, 'NAXIS2')
	naxis = sxpar(h, 'NAXIS')
	if (naxis eq 3) then ns = sxpar(h, 'NAXIS3') else ns = 1
	

; Check to see if BCD or not! (1: BCD, 0: Raw, -1: pbcd?)
	areadmod = sxpar(h, 'AREADMOD', COUNT=acount)
	aortime = sxpar(h, 'AORTIME', COUNT=xcount)
	if (xcount gt 0) then bcdflag = 1 else bcdflag = 0
	if (acount eq 0) then bcdflag = -bcdflag

; Get times
;	if (bcdflag eq 0) then begin
;		areadmod = sxpar(h, 'A0617D00')

; The sclk of the observation is
;		a649 = double(sxpar(h, 'A0649D00'))
;		a650 = double(sxpar(h, 'A0650D00'))
;		a612 = double(sxpar(h, 'A0612D00'))
;		a652 = double(sxpar(h, 'A0652D00'))
;		sclk = a649 + a650 / 65536.D + 0.01D * (a612 - a652)
; The calculation below gives the time that the DCE was received by the S/C
; C+DH
;;		sclk = double(sxpar(h,'H0122D00')) + $
;;		         double(sxpar(h,'H0123D00')) / 256.0
;		arrmode = sxpar(h, 'A0617D00')
;		if (arrmode eq 1) then clock = 0.01 else clock = 0.2
;		ft = clock * (2.0 * sxpar(h, 'A0614D00') + sxpar(h, 'A0615D00'))
;		dt = findgen(ns) * ft
;		t = dt + sclk
;		hjd = t * 0.
;	endif else begin
;		sclk = sxpar(h, 'SCLK_OBS')
		ft = sxpar(h, 'FRAMTIME')
;		dt = (findgen(ns) + 0.5) * ft
;		xft = dblarr(ns) + ft
;		if (areadmod eq 1) then xft = -xft
;		t = dt + sclk
;		hjd0 = sxpar(h, 'HMJD_OBS')
;		hjd = hjd0 + dt / (24.D * 3600.D)
;	endelse

; If raw image need to convert to cube
	if (bcdflag eq 0 and areadmod eq 1) then begin
		nx = 32
		ny = 32
		ns = 64
		
; transform to cube
		cube = reform(im, nx, ny, ns)
; Perform InSb sign flip
		if (nch lt 3) then cube = fix(65535 - cube)
; Rotate to BCD coordinates
		if (nch lt 3) then cube = reverse(cube, 2)
; Blank bad pixel - only know about bad pixel in channel 1 raws so far
		if (nch eq 1 and not keyword_set(WARM)) then cube[24, 6, *] = 0
	endif else cube = im

; Get sigma image
; Need to determine which readnoise and sat level to use
; 0.02, 0.1, 0.4, 2.0 sub, 0.4, 0.6, 1.2, 2, 6, 12, 30, 100 frametimes
	if (abs(ft-0.02) lt 0.01) then findex = 0 $
	else if (abs(ft-0.1) lt 0.01) then findex = 1 $
	else if (abs(ft-0.4) lt 0.01 and areadmod eq 1) then findex = 2 $
	else if (abs(ft-0.4) lt 0.01 and areadmod eq 0) then findex = 4 $
	else if (abs(ft-0.6) lt 0.01) then findex = 5 $
	else if (abs(ft-1.2) lt 0.01) then findex = 6 $
	else if (abs(ft-2.) lt 0.01 and areadmod eq 1) then findex = 3 $
	else if (abs(ft-2.) lt 0.01 and areadmod eq 0) then findex = 7 $
	else if (abs(ft-6.) lt 0.01) then findex = 8 $
	else if (abs(ft-12.) lt 0.01) then findex = 9 $
	else if (abs(ft-30.) lt 0.01) then findex = 10 $
	else if (abs(ft-100.) lt 0.01) then findex = 11 else findex = 0

; If no sigma image then make one -- BCD case first poisson + readnoise
	if (bcdflag ne 0) then begin
		sbtoe = sxpar(h, 'GAIN') * sxpar(h, 'EXPTIME') / sxpar(h,'FLUXCONV')
		sig = abs(cube) * sbtoe
; Make noise positive definite
;			if (min(sig) lt 0.0) then sig = sig - min(sig)
		sig = sig + readnoise[nch-1, findex] * readnoise[nch-1, findex]
		sig = sqrt(sig)
; Convert back to MJy/sr
		sig = sig / sbtoe
	endif else begin
; Sigma image for raw data, poisson noise plus readnoise
		dntoe = gain[nch-1] 
		sig = cube * dntoe
		sig = sig + readnoise[nch-1, findex] * readnoise[nch-1, findex]
		sig = sqrt(sig)
; Convert back to DN
		sig = sig / dntoe
	endelse

; Create arrays to hold centroids
;XX channging doubles to floats; check that this is ok
	x3 = dblarr(ns)* !VALUES.D_NAN
;	y3 = fltarr(ns)* !VALUES.D_NAN
;        x3 = replicate(!VALUES.D_NAN, ns)
        y3 = x3
;	x5 = fltarr(ns) * !VALUES.D_NAN
;	y5 = fltarr(ns)* !VALUES.D_NAN
;	x7 = fltarr(ns)* !VALUES.D_NAN
;	y7 = fltarr(ns)* !VALUES.D_NAN
;	xc = fltarr(ns)* !VALUES.D_NAN
;	yc = fltarr(ns)* !VALUES.D_NAN
;	xg = fltarr(ns)* !VALUES.D_NAN
;	yg = fltarr(ns)* !VALUES.D_NAN
;	xh = fltarr(ns)* !VALUES.D_NAN
;	yh = fltarr(ns)* !VALUES.D_NAN

; Centroids using alternate uncertainty estimate
	xp3 = x3
	yp3 = x3
;	xp5 = x3
;	yp5 = x3
;	xp7 =x3
;	yp7 = x3

; Flux and background arrays, 4 choices of background annulus
	f = dblarr(ns, napers)* !VALUES.D_NAN
	b = f; fltarr(ns, nbacks)* !VALUES.D_NAN
	fp =  f;fltarr(ns, napers)* !VALUES.D_NAN
	bp =  f;fltarr(ns, nbacks)* !VALUES.D_NAN
	bb =  f;fltarr(ns,nbacks) * !VALUES.D_NAN; from box_centroider directly
;  count arrays
	c =  f;fltarr(ns, napers)* !VALUES.D_NAN
	cb =  f;fltarr(ns, nbacks)* !VALUES.D_NAN

; Sigma arrays
	x3s = x3; fltarr(ns)* !VALUES.D_NAN
	y3s = x3;fltarr(ns)* !VALUES.D_NAN
;	x5s = x3;fltarr(ns)* !VALUES.D_NAN
;	y5s = x3;fltarr(ns)* !VALUES.D_NAN
;	x7s = x3;fltarr(ns)* !VALUES.D_NAN
;	y7s = x3;fltarr(ns)* !VALUES.D_NAN

	xp3s = x3;fltarr(ns)* !VALUES.D_NAN
	yp3s = x3;fltarr(ns)* !VALUES.D_NAN
;	xp5s = x3;fltarr(ns)* !VALUES.D_NAN
;	yp5s = x3;fltarr(ns)* !VALUES.D_NAN
;	xp7s = x3;fltarr(ns)* !VALUES.D_NAN
;	yp7s = x3;fltarr(ns)* !VALUES.D_NAN

	fs = f;fltarr(ns, napers)* !VALUES.D_NAN
	bs = f;fltarr(ns, nbacks)* !VALUES.D_NAN
	fps = f;fltarr(ns, napers)* !VALUES.D_NAN
	bps = f;fltarr(ns, nbacks)* !VALUES.D_NAN

; Noise pixel array
	np = dblarr(ns,/nozero)
	
; Subframe array
	sf = lindgen(ns)
	
; Flag array
	flag = intarr(ns)
	
; Full width at half max arrays
	xfwhmarr = dblarr(ns,/nozero)
	yfwhmarr = dblarr(ns,/nozero)

;
; Find centroid for each image plane
	for i = 0, ns-1 do begin
           slice = cube[*, *, i]
;           bslice = bkgd[*,i]
; Uncertainty image calculated from Poisson plus readnoise
           sigma2 = sig[*, *, i] * sig[*, *, i]
; SSC provided uncertainty image		
           unc2 = unc[*, *, i] * unc[*, *, i]

; Find position of brightest pixel or if coordinate is passed find brightest
; pixel in small search window (5 pixel) around coordinate
           skip_src = 0
           adxy, h, ra, dec, xmax, ymax

;cheating for now
;           xmax = xmax
;           ymax = ymax 
; Only continue with images that have a source that not 5 pixels from edge
          ;; xedge = nx - edge - 1.
           xedge = nx - edge + 1.  ;;jk changing to allow closer to edge.
           if (xmax lt edge or xmax gt xedge or ymax lt edge or ymax gt xedge)  then begin ;   or eflux[nch-1] gt sat_levels[nch-1, findex])
		    
		    ;if (eflux[nch-1] gt sat_levels[nch-1, findex]) then begin
		    ;	print, 'Source ' + strn(eflux[nch-1]) + 'brighter than ' + $
		    ;	       strn(sat_levels[nch-1, findex]) + ' for ' + sxpar(h, 'RAWFILE')
		    ;	flag[0:ns-1] = -1
		  ;	endif
              if (xmax lt edge or xmax gt xedge or ymax lt edge or ymax gt xedge) then begin 
                 expid =  sxpar(h, 'EXPID')
                 ;;print, 'xmax, xedge, edge, ymax', xmax, xedge, edge, ymax
                 if not keyword_set(silent) then print, 'Source outside array for image number ', expid
                 flag[0:ns-1] = -2
              endif 

              skip_src = 1
           endif
		
; If source position is in image, then try to find centroid and perform photometry
           if (skip_src eq 0) then begin
; Calculate size of pixels in arcseconds
; Code that used header information so full and subarray frames used different plate scales
; which is wrong as all corrections applied assume a uniform plate scale to account for
; distortions and other variations across the focal plabe
;;			if (bcdflag ne 0) then pscale2 = abs(sxpar(h, 'PXSCAL1')) * $
;;	                                    sxpar(h, 'PXSCAL2') $
;;	        else pscale2 = 1.22D * 1.22D
; Updated code using hardwired constants 11 May 2011 SJC	        
	        pscale2 = pxscal1[nch-1] * pxscal2[nch-1]
	        pscale = sqrt(pscale2)

; Calculate source centroid a variety of ways, will use 3 pixel half radius 1st moment centroider
; as default

; if keyword APER set, then use it for the aperture photometry
; 3 pixel half-width, use smallest aperture for FWHM
;                print, 'starting box_centroider ', xmax, ymax
                expid = sxpar(h, 'EXPID')
                box_centroider, slice, sigma2, xmax, ymax, 3, 6,3, tx, ty, tf, tb, $
                                txs, tys, tfs, tbs, tc, tcb, tnp, xfwhm, yfwhm, expid,/twopass

;pro box_centroider, input_image, sigma2, xmax, ymax, halfboxwidth, $
;                    backboxwidth, boxborder, x0, y0, f0, b, xs, ys, fs, bs, $
;                    c, cb, np, xfwhm, yfwhm, xys, xycov, expid, MMM=mmm



;                print, 'finished box_centroider', tx, ty
                x3[i] = tx & y3[i] = ty & x3s[i] = txs & y3s[i] = tys
                bb[i] = tb
                np[i] =tnp
; Store FWHM			
                xfwhmarr[i] = xfwhm & yfwhmarr[i] = yfwhm
; 5 pixel half-width
;			box_centroider, slice, sigma2, xmax, ymax, 5, 6, 3, tx, ty, tf, tb, $
;			                txs, tys, tfs, tbs, tc, tcb, tnp
;			x5[i] = tx & y5[i] = ty & x5s[i] = txs & y5s[i] = tys
; 7 pixel half-width
;			box_centroider, slice, sigma2, xmax, ymax, 7, 6, 3, tx, ty, tf, tb, $
;			                txs, tys, tfs, tbs, tc, tcb, tnp
;			x7[i] = tx & y7[i] = ty & x7s[i] = txs & y7s[i] = tys

; Do same calculation with BCD uncertainty images
; 3 pixel half-width
;                box_centroider, slice, unc2, xmax, ymax, 3, 6, 3, tx, ty, tf, tb, $
;                                txs, tys, tfs, tbs, tc, tcb, tnp
;                xp3[i] = tx & yp3[i] = ty & xp3s[i] = txs & yp3s[i] = tys
; 5 pixel half-width
;			box_centroider, slice, unc2, xmax, ymax, 5, 6, 3, tx, ty, tf, tb, $
;			                txs, tys, tfs, tbs, tc, tcb, tnp
;			xp5[i] = tx & yp5[i] = ty & xp5s[i] = txs & yp5s[i] = tys
; 7 pixel half-width
;			box_centroider, slice, unc2, xmax, ymax, 7, 6, 3, tx, ty, tf, tb, $
;			                txs, tys, tfs, tbs, tc, tcb, tnp
;			xp7[i] = tx & yp7[i] = ty & xp7s[i] = txs & yp7s[i] = tys
			
; Use IDLphot centroid methods
; Gaussian centroider
;			gcntrd, slice, xmax, ymax, tx, ty, fwhm[nch-1] / pscale, SILENT=silent
;			xg[i] = tx & yg[i] = ty
; Other centroider
;			cntrd, slice, xmax, ymax, tx, ty, fwhm[nch-1] / pscale, SILENT=silent
;			xh[i] = tx & yh[i] = ty
		
; Convert image to electrons
                if (bcdflag ne 0) then begin
                   eim = slice * sbtoe 
                endif else begin
                   eim = slice
                endelse


;manually find the background in the images
;mask the top and bottom row
;                eim[*,0] = !VALUES.D_NAN &  eim[*,31] = !VALUES.D_NAN
;                ebslice = calc_bkgd (eim, h, ra, dec)

; 1st set of apertures
;;			aper, eim, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps1, back1, $
;;                aper, eim, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps1, back1, $
                aper, eim, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps1, back6, $
                      badpix, /FLUX, /EXACT, /NAN, /SILENT, /MEANBACK,$
                      READNOISE=readnoise[nch-1, findex];, SETSKYVAL = ebslice
                f[i, 0:(naps1-1)] = xf / sbtoe
                b[i, 0] = xb
                bs[i, 0] = xbs
                bptr = where(xf ne xf, bcount)
                if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
                rn = readnoise[nch-1, findex] * !DPI * aps1
                fs[i, 0:(naps1-1)] = sqrt(xfs * xfs + rn * rn) / sbtoe

; 2nd set of apertures
;			aper, eim, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps2, back2, $
;                aper, eim, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps2, back4, $  ;
;                      badpix, /FLUX, /EXACT, /SILENT, /NAN, /meanback,$         ;
;                      READNOISE=readnoise[nch-1, findex]                        ;
;                f[i, naps1] = xf / sbtoe                                        ;
;                b[i, 1] = xb
;                bs[i, 1] = xbs
;                bptr = where(xf ne xf, bcount)
;                if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;                rn = readnoise[nch-1, findex] * !DPI * aps2
;                fs[i, naps1] = sqrt(xfs * xfs + rn * rn) / sbtoe

; 3rd set of apertures
;			aper, eim, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps3, back3, $
;                aper, eim, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps2, back5, $
;                      badpix, /FLUX, /EXACT, /SILENT, /NAN, /MEANBACK,$
;                      READNOISE=readnoise[nch-1, findex]
;                                ;                      ;print, 'xf', xf
;                f[i, naps1+1] = xf / sbtoe
;                b[i, 2] = xb
;                bs[i, 2] = xbs
;                bptr = where(xf ne xf, bcount)
;                if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;                rn = readnoise[nch-1, findex] * !DPI * aps2
;                fs[i, naps1+1] = sqrt(xfs * xfs + rn * rn) / sbtoe

; 4th set of apertures
;			aper, eim, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps4, back4, $
;			aper, eim, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps4, back4, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN, /MEANBACK, $
;			      READNOISE=readnoise[nch-1, findex]
;			f[i, 10] = xf / sbtoe
;			b[i, 3] = xb
;			bs[i, 3] = xbs
;			bptr = where(xf ne xf, bcount)
;			if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;			rn = readnoise[nch-1, findex] * !DPI * aps4
;			fs[i, 10] = sqrt(xfs * xfs + rn * rn) / sbtoe

; Now using BCD uncertainty image to determine noise by summing uncertainty^2 in same
; aperture, in this case uncertainty is the uncertainty summed over the aperture plus
; the uncertainty in the background subtracted.

; 1st set of apertures
;;			aper, slice, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps1, back1, $
;                aper, slice, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps1, back1, $
;                      badpix, /FLUX, /EXACT, /SILENT, /NAN, /MEANBACK
;                fp[i, 0:(naps1 - 1)] = xf
;                bptr = where(xf ne xf, bcount)
;;			aper, unc2, x5[i], y5[i], xfs, txfs, txb, txbs, 1.0, aps1, back1, $
;                aper, unc2, x3[i], y3[i], xfs, txfs, txb, txbs, 1.0, aps1, back1, $
;                      badpix, /FLUX, /EXACT, /SILENT, /NAN,/MEANBACK
;                if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;                nsource = !DPI * aps1 * aps1
;                nback = !DPI * (back1[1] * back1[1] - back1[0] * back1[0])
;                fps[i, 0:(naps1-1)] = xfs + xbs * nsource * nsource / nback
                
; 2nd set of apertures
;;			aper, slice, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps2, back2, $
;			aper, slice, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps2, back2, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN,/MEANBACK
;			fp[i, 8] = xf
;			bptr = where(xf ne xf, bcount)
;;			aper, unc2, x5[i], y5[i], xfs, txfs, txb, txbs, 1.0, aps2, back2, $
;			aper, unc2, x3[i], y3[i], xfs, txfs, txb, txbs, 1.0, aps2, back2, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN,/MEANBACK
;			if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;			nsource = !DPI * aps2 * aps2
;			nback = !DPI * (back2[1] * back2[1] - back2[0] * back2[0])			
;			fps[i, 8] = xfs + xbs * nsource * nsource / nback

; 3rd set of apertures
;;			aper, slice, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps3, back3, $
;			aper, slice, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps3, back3, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN,/MEANBACK
;			fp[i, 9] = xf
;			bptr = where(xf ne xf, bcount)
;;			aper, unc2, x5[i], y5[i], xfs, txfs, txb, txbs, 1.0, aps3, back3, $
;			aper, unc2, x3[i], y3[i], xfs, txfs, txb, txbs, 1.0, aps3, back3, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN, /MEANBACK
;			if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;			nsource = !DPI * aps3 * aps3
;			nback = !DPI * (back3[1] * back3[1] - back3[0] * back3[0])					
;			fps[i, 9] = xfs + xbs * nsource * nsource / nback

; 4th set of apertures
;			aper, slice, x5[i], y5[i], xf, xfs, xb, xbs, 1.0, aps4, back4, $
;			aper, slice, x3[i], y3[i], xf, xfs, xb, xbs, 1.0, aps4, back4, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN
;			fp[i, 10] = xf
;			bptr = where(xf ne xf, bcount)
;;			aper, unc2, x5[i], y5[i], xfs, txfs, txb, txbs, 1.0, aps4, back4, $
;			aper, unc2, x3[i], y3[i], xfs, txfs, txb, txbs, 1.0, aps4, back4, $
;			      badpix, /FLUX, /EXACT, /SILENT, /NAN
;			if (bcount ne 0) then xfs[bptr] = !VALUES.D_NAN
;			nsource = !DPI * aps4 * aps4
;			nback = !DPI * (back4[1] * back4[1] - back4[0] * back4[0])						
;			fps[i, 10] = xfs + xbs * nsource * nsource / nback
; Set flag to number of good apertures
			ptr = where(f[i, *] eq f[i, *], count)
			flag[i] = count					
		endif else pscale2 = 1.0
	endfor
	
; Convert from summed MJy/sr to Jy
; convert scale from arcsec^2 to sr and scale to Jy
	scale = pscale2 * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06
	f = f * scale
	fp = fp * scale
	b = b * scale
        bb = b * scale
; Return sigma instead of variance
	x3s = sqrt(x3s)
	y3s = sqrt(y3s)
;	x5s = sqrt(x5s)
;	y5s = sqrt(y5s)
;	x7s = sqrt(x7s)
;	y7s = sqrt(y7s)

	xp3s = sqrt(xp3s)
	yp3s = sqrt(yp3s)
;	xp5s = sqrt(xp5s)
;	yp5s = sqrt(yp5s)
;	xp7s = sqrt(xp7s)
;	yp7s = sqrt(yp7s)

	fs = fs * scale
	fps = sqrt(fps)
	fps = fps * scale
	bs = bs * scale
	
; If subarray mode, return the equivalent full array pixel coordinates
;	if (areadmod eq 1) then begin
;		rowoff = 8.0D
;		x3 = x3 + rowoff
;		x5 = x5 + rowoff
;		x7 = x7 + rowoff
;		xg = xg + rowoff
;		xh = xh + rowoff
;		xp3 = xp3 + rowoff
;		xp5 = xp5 + rowoff
;		xp7 = xp7 + rowoff	
		
;		if (nch lt 3) then coloff = 216.D else coloff = 8.D
;		y3 = y3 + coloff
;		y5 = y5 + coloff
;		y7 = y7 + coloff
;		yg = yg + coloff
;		yh = yh + coloff
;		yp3 = yp3 + coloff
;		yp5 = yp5 + coloff
;		yp7 = yp7 + coloff	
;	endif

return
end


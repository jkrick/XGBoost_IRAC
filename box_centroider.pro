pro box_centroider, input_image, sigma2, xmax, ymax, halfboxwidth, $
                    backboxwidth, boxborder, x0, y0, f0, b, xs, ys, fs, bs, $
                    c, cb, np, xfwhm, yfwhm, xys, xycov, TWOPASS=twopass, MMM=mmm, verbose = verbose
;+
; NAME:
;    BOX_CENTROIDER
; 
; PURPOSE:
;    Calculates centroids for a source using a simple 1st moment box centroider
;    Designed to work with IRAC bcds and uncertainty images, it should work
;    with other data but has not been tested with such.
; 
; INPUTS:
;    INPUT_IMAGE: 2d float/double array containing input image, input image is 
;           presumed to be have bad pixels NaNed out       
;    SIGMA2: 2d float/double array containing estimate of square of per pixel  
;            image uncertainty
;    XMAX: float/double scalar containing x position of peak pixel to 
;          centroid around
;    YMAX: float/double scalar containing y position of peak pixel to 
;           centroid around
;    HALFBOXWIDTH: integer scalar, for size of centroiding box, 
;              actually halfboxwidth - 0.5 as aperture size is 2*HALFBOXWIDTH+1
;              Odd number of pixels is used so aperture is centered on 
;              peak pixel
;    BACKBOXWIDTH: integer scalar, for size of rectangular background box
;    BOXBORDER: integer scalar, offset of background annulus in pixels from 
;               box aperture
;
; CONTROL KEYWORD:
;    /MMM: if set, then backround estimate will be calculated using mmm.pro
;    /TWOPASS: if set run this centroider once to find the pixel that holds 
;              the centroid, then rerun using that pixel.
;
; OUTPUTS:
;    X0: double scalar X position of centroid
;    Y0: double scalar Y position of centroid
;    F0: double scalar flux of source as summed in box after background removal
;    B: double scalar background level
;    XS: double scalar uncertainty in X centroid position
;    YS: double scalar uncertainty in Y centroid position
;    FS: double scalar uncertainty in measured flux accounting for uncertainty 
;        in background
;    BS: double scalar uncertainty in background
;    C: fixed scalar number of good pixels in aperture
;    CB: fixed scalar number of good pixels in background aperture
;    NP: fixed scalar number of noise pixels
;    XFWHM: full width at half maximum in x-coordinate of centroid
;    YFWHM: full width at half maximum of y-coordinate of centroid
;    XYS: double scalar covariance in X wrt Y centroid position (and vice 
;            versa) (Propagated uncertainty)
;    XYCOV: x-y covariance of centroid (statistics on input data) 
;
;
; ALGORITHM:
;    The centroid position in units of pixels on the image is calculated
;    using flux-weighted 1st moments in x and y pixel position.  
;          x0 = sum(I(i) * x(i)) / sum(I(i)) 
;          y0 = sum(I(i) * y(i)) / sun(I(i))
;    where the sum is conducted over all pixels i in the source aperture, 
;    a square box of size 2*HALFBOXWIDTH+1 pixels centered on pixel (XMAX,YMAX).
;
;    Before calculating the moments, an estimate of the background is determined ;    from a rectangular background region which is a specified width 
;    (BACKBOXWIDTH) and offset (BOXBORDER) from the source aperture used to 
;    calculate the moments.  By default, a sigma-clipped median is used for the
;    background value.
;
;    The uncertainties and covariance for the centroid position are calculated 
;    from the supplied per-pixel uncertainty by propagating the per-pixel 
;    uncertainties:
;		xs^2 = sum((x(i)-x0)^2 * sigma(i)^2)
;		ys^2 = sum((x(i)-x0)^2 * sigma(i)^2)
;       cov = sum((x(i)-x0) * (y(i)-y0) * sigma(i)^2)
;    where the sums are again over all pixels i in the source aperture
;
;    Source FWHM is calculated from the 2nd moment of light:
;       xfwhm = 2.D * sqrt(2.D * alog(2.D)) * 
;                         [sum(x(i)-x0)^2 * I(i)) / sum(I(i))]^0.5
;       yfwhm = 2.D * sqrt(2.D * alog(2.D)) * 
;                         [sum(y(i)-y0)^2 * I(i)) / sum(I(i))]^0.5
;
;     and a measure of the size of the source in terms of noise pixels is
;     calculated:
;       noise pixels = sum(I(i))^2 / sum(I(i)^2)
;     where the default aperture size should be sufficient for IRAC data.
;     
; NOTES:
;    The aperture size and background regions have been optimized for IRAC
;    BCDs and should be modified accordingly for other data. 

;    1st moment centroiding provides a robust measurement of source position
;    for undersampled data where more complex models fail due to lack of
;    sampling or inadequate knowledge of the point response function.
;
;    Any centroiding method does have inherent biases which should be 
;    understood when interpreting the positions, particularly if the measured 
;    flux depends on the location of the source relative to a pixel center
;    due to the effects of intra-pixel gain variations.
;
;    The flux returned is a more crude estimate than achieved with
;    circular/elliptical aperture photometry that uses fractional pixels and is
;    only provided for reference and sanity checking.
;
; HISTORY:
;    Added correct np output for the case where 50% of pixels in source aperture are bad  10 jul 14 JK
;    Changed shape of background aperture to be square 6 May JK
;    Added keyword /TWOPASS to find the pixel that holds the centroid - 
;        avoids discrepancies based on where the box is initially centered 
;        6 May 2014 JGI
;    Use ROUND() to set the box limits 17 Apr 2014 JGI
;    Set f0 = fluxsum for return  13 Aug 2013 JGI 
;    Added additional header documentation. 09 Aug 2013 SJC
;    Added minimum number of pixels check to calculate 
;    background 09 Aug 2013 SJC
;    Added outputs XYS, XYCOV   14 Jun 2013  JGI
;    Corrected calculation of noise pixels 21 June 2013 SJC
;    Need to not modify the image directly by subtracting the background as 
;    that corrupts the background level 12 Feb 2013 SJC
;    Added source width calculation 27 July 2011 SJC
;    Recoding from get_centroids_nf procedure 26 May 2010 SJC
;    
;-
  IF KEYWORD_SET(TWOPASS) THEN BEGIN
     ;;print, 'twopass'
     box_centroider, input_image, sigma2, xmax, ymax, halfboxwidth, $
                     backboxwidth, boxborder, newxmax, newymax, MMM=mmm
     if newxmax gt 0 and newymax gt 0 then begin
        xmax = newxmax
        ymax = newymax
     endif
     ;;if keyword_set(verbose) print, 'xmax, ymax ', xmax, ymax
 
  ENDIF
; Copy input image to working image
  image = input_image
 
; Set sigma clipping threshold
	sigma = 3

; Get image size
	sz = size(image)
	nx = sz[1]
	ny = sz[2]
		
; Arrays holding x and y coordinates of each pixel -- allow arbitrary sized 
; images
	xx = findgen(nx) # replicate(1.0, ny)
	yy = replicate(1.0, nx) # findgen(ny)
	
; Set box limits for aperture
	xa = ROUND(xmax-halfboxwidth) > 0 < (nx-1)
	xb = ROUND(xmax+halfboxwidth) < (nx-1) > 0
	ya = ROUND(ymax-halfboxwidth) > 0 < (ny-1)
	yb = ROUND(ymax+halfboxwidth) < (ny-1) > 0
	npixels_in_box = (2*halfboxwidth+1) * (2*halfboxwidth+1)	

; Set background regions
	xamin = xa - backboxwidth - boxborder > 0 < xa
	xamax = xa - boxborder > 0 < xa
	xbmax = xb + backboxwidth + boxborder > xb < (nx-1)
	xbmin = xb + boxborder > xb < (nx-1)
	yamin = ya - backboxwidth - boxborder > 0 < ya
	yamax = ya - boxborder > 0 < ya
	ybmax = yb + backboxwidth + boxborder > yb < (ny-1)
	ybmin = yb + boxborder > yb < (ny-1)

; Segment image into box (mask = 1), background (mask=-1) and 
; everything else (mask = 0)
	mask = image * 0
	mask[xa:xb, ya:yb] = 1
	bptr = where(xx le xamax and xx ge xamin and yy le ybmax and yy ge yamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(xx le xbmax and xx ge xbmin and yy le ybmax and yy ge yamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(yy le yamax and yy ge yamin and xx le xbmax and xx ge xamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(yy le ybmax and yy ge ybmin and xx le xbmax and xx ge xamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
; 


; remove background - median of background region, unless MMM is used
	bptr = where(mask eq -1 and image eq image, bcount)
	min_good_background_pixels = 0.5 * npixels_in_box
	if (bcount ge min_good_background_pixels) then begin 
		if (keyword_set(MMM)) then begin
			mmm, image[bptr], back, bsig
			bsig = bsig * bsig
		endif else begin
			back = median(image[bptr]) 
; One pass of outlier rejection for the background, could use mmm as well
			bmom = moment(image[bptr])
			bsig = bmom[1]
			bbptr = where(abs(image[bptr] - back) le sigma * sqrt(bsig), $
					               bbcount)
			if (bbcount ge min_good_background_pixels) then begin
				back = median(image[bptr[bbptr]])
				bmom = moment(image[bptr[bbptr]] - back)
				bsig = bmom[1]
				bcount = bbcount
			endif else begin
				back = 0.0D
				bsig = 0.0D
			endelse
		endelse
	endif else begin
		back = 0.0D
		bsig = 0.0D
	endelse

; Perform background subtraction
	image = image - back
	b = back

	min_good_pixels_in_aperture = 0.5 * npixels_in_box
; Perform centroid in 2*boxwidth+1 pixel box centered on peak pixel
        gptr = where(mask eq 1 and image eq image, gcount)
;;        print,'gptr', image[gptr]
; Only calculate centroid if more than 50% of pixels in source aperture are good	
	if (gcount gt min_good_pixels_in_aperture) then begin
		fluxsum = total(image[gptr])
		flux2sum = total(image[gptr] * image[gptr])
		fluxsum2 = fluxsum * fluxsum
		x0 = total(xx[gptr] * float(image[gptr])) / fluxsum
                y0 = total(yy[gptr] * float(image[gptr])) / fluxsum
		f0 = fluxsum   ;; Added 13 Aug 2013 JGI
; Now calculate the variance 
; note for normal distributions, fwhm = 2*sqrt(2 * ln 2) sigma ~ 2.35482 sigma
		dx = xx[gptr] - x0
		dy = yy[gptr] - y0
		xfwhm = total(dx * dx * float(image[gptr])) / fluxsum
		yfwhm = total(dy * dy * float(image[gptr])) / fluxsum
; Added covariance 09 Aug 2013 JGI
		xycov = total(dx * dy * float(image[gptr])) / fluxsum
		sigscale = 2.D * sqrt(2.D * alog(2.D))
		xfwhm = sigscale * sqrt(xfwhm)
		yfwhm = sigscale * sqrt(yfwhm)
		
; Use summed flux as measure of total flux and calculate noise pixels
; Need a measure of the PRF as determined from the image
; so we can use I(x,y) = F * P(x,y) where I(x,y) is the per pixel flux density,
; F is the flux density and P(x,y) is the PRF value for pixel x,y
; Number of noise pixels is N = 1. / sum(P^2) but P = I/F so 
; N = 1. / sum(I^2/F^2)	= F^2 / sum(I^2) = fluxsum2 / flux2sum
; Updated 21 June 2013 SJC
		np = fluxsum2 / flux2sum

; Number of pixels in source and background regions
		c = float(gcount)
		cb = float(bcount)
		
; If a bad pixel exists in the source aperture return a flux of NaN.
		xptr = where(mask eq 1 and image ne image, xcount)
		if (xcount gt 0) then f = !VALUES.D_NAN else f = fluxsum
	
; Uncertainties
		xs = total(dx * dx * sigma2[gptr])
		ys = total(dy * dy * sigma2[gptr])
; Covariance
		xys = total(dx * dy * sigma2[gptr])
		xs = xs / fluxsum2
		ys = ys / fluxsum2
		xys = xys / fluxsum2
		bs = c * c * bsig / cb
		fs = total(sigma2[gptr]) + bs
	endif else begin
		print, 'BOX_CENTROIDER: Less than 50% good pixels in source box'
		x0 = !VALUES.D_NAN
		y0 = !VALUES.D_NAN
		f0 = !VALUES.D_NAN  ;; Updated 13 Aug 2013 JGI
		b = !VALUES.D_NAN
		xs = !VALUES.D_NAN
		ys = !VALUES.D_NAN
		xys = !VALUES.D_NAN
		bs = !VALUES.D_NAN
		fs = !VALUES.D_NAN
		xfwhm = !VALUES.D_NAN
		yfwhm = !VALUES.D_NAN
		xycov = !VALUES.D_NAN
                np = !VALUES.D_NAN
	endelse

return
end

pro dewrap2,im,ichan,barrel,fowlnum,pedsig,pattern,im_out
; 2003-03-29 BG takes raw image and removes wraparound and converts sign.
; makes best estimate of the truncation error and corrects for
; it where appropriate. Assumes for Fowler frames that the most
; negative true fowler signal (after the barrel shift) is -8000 DN.
; NOTE: Remember that when computing a Fowler-equivalent, which is
;   not done in this routine, for channels 1 and 2, it's P - S and
;   for channels 3 and 4 it's S - P.
;
; for all standard frames, fowlnum = 2^barrel -- that is, the barrel
; shift turns the accumulation into an average
;
;INPUT
; im (fixed: read in as 16 bit signed integers)
; ichan	IRAC channel number 1-4
; barrel ABARREL barrel shift
; fowlnum AFOWLNUM
; pedsig APEDSIG 1=pedesal 2= signal 0 otherwise
; pattern 0=no pattern  1-4 = patterns 1-4
;output  
; im_out  (floating)
;
;**********************************************************
min_possible_value= -8000. ;lowest value a Fowler frame is allowed to have --
;                 lower values will be assumed to be higher than
;                 32767. The maximum possible value is 
;		  65535 +  min_possible_value
;		  This is done after correcting for truncation, barrel
;		  shift, and fowler number, and rectification of channels
;		  1 and 2.
;
im_out=float(im)
if fowlnum eq 0 then begin      ;FOWLER-0 and PATTERNs are in range 0-65535
;		patterns are ramp 0 to FFFFH and ramp FFFFH to 0 and
;               checkerboards 5555H (pos) and AAAAH (neg)
;		Fowler-0 or fixed patterns require only the correction for
;			signed integer to unsigned integer    
   a=where(im lt 0)
   if a[0] ne -1 then im_out[a]=im_out[a]+65536.
   return
endif
if pedsig eq 1 then begin       ;PEDESTAL FRAME
  a=where(im gt 0)
  if a[0] ne -1 then im_out[a]=im_out[a]-65536.
;  correct for bit truncation 
  corr1=0.5*(1.-2.^(-barrel))
  if barrel ne 0 then im_out=im_out+corr1
;  correct for barrel shift and fowler multiplication
  if fowlnum ne 2^barrel then im_out=im_out*(2.^barrel/fowlnum)
;   negate the pedestal -- it was loaded negative in DSP
  im_out= -im_out	;negate pedestal since it was loaded negative  
  return
endif else if pedsig eq 2 then begin    ;SIGNAL FRAME
  a=where(im lt 0)
  if a[0] ne -1 then im_out[a]=im_out[a]+65536.
;  correct for bit truncation 
  corr1=0.5*(1.-2.^(-barrel))
  if barrel ne 0 then im_out=im_out+corr1
;  correct for barrel shift and fowler multiplication
  if fowlnum ne 2^barrel then im_out=im_out*(2.^barrel/fowlnum)
  return
endif else begin                ;FOWLER FRAME
;  correct for bit truncation 
  corr1=0.5*(1.-2.^(-barrel))
  if barrel ne 0 then im_out=im_out+corr1
;  correct for barrel shift and fowler multiplication
  if fowlnum ne 2^barrel then im_out=im_out*(2.^barrel/fowlnum)
; negate the values for channels 1 and 2
  if ichan eq 1 or ichan eq 2 then im_out= -im_out
; now remove the assumed wraparound
  a=where(im_out lt min_possible_value)
  if a[0] ne -1 then im_out[a]=im_out[a]+65536.
  return
endelse
end

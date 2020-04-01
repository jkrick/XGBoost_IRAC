 ; cernox temperature conversion
; p.131 of IRAC-582-SPEC-004
;
; RANGE=1 is the default CRYO
; RANGE=4 is patched POST-CRYO
;
function fpa1t2,dn,$
                range,$ ; cernox range, required but ignored if OHMS
                OHMS=OHMS,$ ; input dn are actually ohms
                DIODE=DIODE,$ ; return a diode temp instead of cernox
                VERBOSE=VERBOSE ; print some stuff


coeff=fltarr(8)
coeff[*]=0

dn=dn*1.0 ; promote to float

if ( keyword_set(OHMS) NE 1) then begin

r = 2.842632998e-2 * (dn+32767) + 4.717189028


case 1 of
(range EQ 1): factor=1.0
(range eq 2): factor= 40.031/107.406
(range eq 3): factor = 40.031/15.677
(range eq 4): factor= 40.031/85.00
else: factor=0
endcase


r=r*factor

endif else begin

r=dn

endelse

if keyword_set(VERBOSE) then print,'R (ohms): ',r


case 1 of
; 16.01 - 88.98 K
((r LE 689) AND (r GT 143)): begin 
   coeff[0]= 3.801895464e2
   coeff[1]= -4.431111640
   coeff[2]= 2.702000353e-2
   coeff[3]= -9.748935274e-5
   coeff[4]= 2.154993624e-7
   coeff[5]= -2.866140516e-10
   coeff[6]= 2.104924192e-13
   coeff[7]= -6.555009122e-17
   outofrange = 0
  end
; 6.004 - 15.983 K
((r LE 2037) AND (r GT 689)): begin
   coeff[0]= 1.123198839e2
   coeff[1]= -4.293536773e-1
   coeff[2]= 8.538518095e-4
   coeff[3]= -9.938826032e-7
   coeff[4]= 7.043196334e-10
   coeff[5]= -2.992786047e-13
   coeff[6]= 7.011393062e-17
   coeff[7]= -6.963814125e-21
   outofrange = 0
  end
; 2.644 - 5.996 K
((r LE 9795) AND (r GT 2037)): begin
   coeff[0]= 2.333640639e1
   coeff[1]= -2.019731650e-2
   coeff[2]= 9.665507982e-6
   coeff[3]= -2.669571437e-9
   coeff[4]= 4.341263936e-13
   coeff[5]= -4.4361858613e-17
   coeff[6]= 2.344914907e-21
   coeff[7]= -5.300647913e-26
   outofrange = 0
  end
endcase


case 1 of
(outofrange EQ 1): tout=9999
(outofrange EQ 0): begin
        x=r
        tout=poly(x,coeff)
     end
else: tout=-1
endcase

if (keyword_set(DIODE) ne 1) then begin
         return,tout
    endif else begin
    
    ; this cernox to diode conversion comes from the 
    ;   v500 and temp setpoint test
    tdiode=[25.38,27.45,29.41,31.36,33.29,35.2,37.1,38.98,40.84]
    tcernox=[24.7,26.2,27.8,29.5,31.4,33.2,35.1,37.1,39]

    
    if (tout LE 25) then return,tout
    if (tout GT 25) then begin

    
    return,interpol(tcernox,tdiode,tout)
    
    
    endif
    
    endelse

end
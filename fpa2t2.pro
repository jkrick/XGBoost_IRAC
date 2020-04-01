; cernox temperature conversion
; p.131 of IRAC-582-SPEC-004
;
; RANGE=1 is the default CRYO
; RANGE=4 is patched POST-CRYO
;
function fpa2t2,dn,$
                range,$ ; cernox range, required but ignored if OHMS
                OHMS=OHMS,$ ; input dn are actually resistance in ohms
                DIODE=DIODE,$ ; return a diode temp instead of cernox
                VERBOSE=VERBOSE ; print some stuff
                
                
coeff=fltarr(8)
coeff[*]=0

if ( keyword_set(OHMS) NE 1) then begin

dn=dn*1.0 ; promote to float

r = 2.650499328e-2 * (dn+32767) + 4.121221611


case 1 of
(range EQ 1): factor=1.0
(range eq 2): factor = 40.037/110.581
(range eq 3): factor = 40.037/15.084
(range eq 4): factor = 40.037/68.
else: factor=0
endcase


r=r*factor

endif else begin

r=dn

endelse

if keyword_set(VERBOSE) then print,'R (ohms): ',r

case 1 of
; 17.037 - 63.818 K
((r LE 646) AND (r GT 248)): begin 
   coeff[0]= 7.554886717e2
   coeff[1]= -9.261701469
   coeff[2]= 5.571612621e-2
   coeff[3]= -1.954933013e-4
   coeff[4]= 4.189419393e-7
   coeff[5]= -5.405833169e-10
   coeff[6]= 3.861064556e-13
   coeff[7]= -1.173214937e-16
   outofrange = 0
  end
; 6 - 16.97 K
((r LE 1438) AND (r GT 646)): begin
   coeff[0]= 4.327501320e2
   coeff[1]= -2.638312089
   coeff[2]= 7.370492818e-3
   coeff[3]= -1.160818400e-5
   coeff[4]= 1.096530024e-8
   coeff[5]= -6.176438515e-12
   coeff[6]= 1.916126254e-15
   coeff[7]= -2.523122090e-19
   outofrange = 0
  end
; 1.849 - 9.938 K
((r LE 5119) AND (r GT 1438)): begin
   coeff[0]= 4.826801515e1
   coeff[1]= -8.600515934e-2
   coeff[2]= 7.552443500e-5
   coeff[3]= -3.792449064e-8
   coeff[4]= 1.144910420e-11
   coeff[5]= -2.056032360e-15
   coeff[6]= 2.024242583e-19
   coeff[7]= -8.412295843e-24
   outofrange = 0
end
else: outofrange = 1
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
    tdiode=[23.4,28.82,30.73,32.64,34.54,36.42,38.3,40.18]
    tcernox=[23.1,27.2,28.7,30.3,31.8,33.3,35.0,36.5]
    
    if (tout LE 25) then return,tout
    if (tout GT 25) then begin

    
    
    return,interpol(tcernox,tdiode,tout)
    
    
    endif
    
    endelse

end

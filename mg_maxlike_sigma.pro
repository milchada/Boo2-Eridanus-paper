;--------------------------------------------------------------------
;+
; NAME:
;     mg_maxlike_sigma
; 
; PURPOSE:
;     Determine central velocity and velocity dispersion based
;     on a maximum likelihood analysis which includes Gaussian
;     distribution on errors on each point.
;
;     Based on expressions from Walker et al (2006, AJ 131, 2114)     
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;     velocity, velocity error
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; MODIFICATION HISTORY:
;   mg  12/06
;- 
;--------------------------------------------------------------------
Pro mg_maxlike_sigma, vobs, sobs, vmax, smax,err_v,err_sig

  
  N = n_elements(vobs)

  ; Need to maximize over true velocity and true sigma
  nvel = 100
  nsig = 100
  lnP   = fltarr(nvel,nsig)

  ; SET UP INITIAL RANGE OVER WHICH TO SEARCH
  vavg = avg(vobs)
  vreal = vavg +  (40 * findgen(nvel)/nvel - 20.0D)
  sreal = 20.0D * findgen(nsig)/nsig 
  Pmax = -1.e6
  for i=0,nvel-1 do begin
     for j=0,nsig-1 do begin
  

       ; USE EQN 8, Walker et al (2006)
        term1  = -0.5D * total( alog(sobs^2 + sreal(j)^2))

        term2a = (vobs - vreal(i))^2
        term2b = (sobs^2 + sreal(j)^2)
        term2  = -0.5D * total(term2a/term2b)

        term3  = -0.5D * N * alog(2/!pi)

        lnP(i,j) = term1 + term2 + term3  


        if lnP(i,j) gt Pmax then begin
           smax = sreal(j)
           vmax = vreal(i)
           Pmax = lnP(i,j)
        endif

     endfor
  endfor   



  ; REFINE SEARCH
  nsig=500
  lnP   = fltarr(nvel,nsig)
  vavg = avg(vobs)
  vreal = vmax +  (5. * findgen(nvel)/nvel - 2.5D)
  sreal = smax * (1. + 0.9*(findgen(nsig)/nsig  - 0.5))
  Pmax = -1e5

  for i=0,nvel-1 do begin
     for j=0,nsig-1 do begin
  

       ; USE EQN 8, Walker et al (2006)
        term1  = -0.5D * total( alog(sobs^2 + sreal(j)^2))

        term2a = (vobs - vreal(i))^2
        term2b = (sobs^2 + sreal(j)^2)
        term2  = -0.5D * total(term2a/term2b)

        term3  = -0.5D * N * alog(2/!pi)

        lnP(i,j) = term1 + term2 + term3  


        if lnP(i,j) gt Pmax then begin
           smax = sreal(j)
           vmax = vreal(i)
           Pmax = lnP(i,j)
        endif

     endfor
  endfor   




  d2lnpu = 0.5 * total(2./(sobs^2 + smax^2))


  o2 = sobs^2 + smax^2
  d2lnpv = -0.5 * total((2./o2) - (4.*smax^2/o2^2)) + $
            0.5 * total((vobs-vmax)^2 * ( (2./o2^2) - (8.*smax^2/o2^3)))

  d2duds = 2.0 * total( (smax*(vobs-vmax)/o2^2))


  detA = (d2lnpu*d2lnpv - d2duds^2)


print,'ML method = '
  err_v  = sqrt(d2lnpv/detA)
  print,'V   = ', vmax,' +/- ',err_v

  err_sig = sqrt(-1.*d2lnpu/detA)
  print,'Sig = ',smax,' +/- ',err_sig


end

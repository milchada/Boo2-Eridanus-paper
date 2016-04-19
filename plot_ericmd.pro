Pro plot_ericmd,create_ps=create_ps

; read MW
 readcol,'data/mw_eri.dat',d,mv,c,t,lte,lg,a,m,ur,mw_gr,ri,iz,mw_r,mux,muy,vr
 plothist,vr,mw_a,mw_b,bin=5


loadct,12
!p.multi=[0,2,1]
!p.charsize=1.1
!p.charthick=1.5

if keyword_set(create_ps) then begin
   set_plot,'ps'
   device,xoffset=2,yoffset=10,filename='fig_ericmd.eps',/times, $
    /portrait, xsize=20,ysize=10,/encapsulate,/color
endif


  ; READ GAL INFO
    maskname = 'Eri'
    path = '~/Projects/Objects/'
    ginfo = mrdfits(path+'pro/galinfo.fits',1)
    q=where(strcompress(ginfo.galname,/rem) eq strcompress(maskname,/rem))
  ; GENERAL GALAXY INFO
    rcen =  ginfo[q].ra
    dcen =  ginfo[q].dec
    dist =  ginfo[q].dist
    EBV =   ginfo[q].EBV


    dmod=5.*alog10(dist*1.e3)-5.
    Ag = 3.793*EBV
    Ar = 2.751 * EBV




 ; READ DEIMOS DATA
    file = 'alldata_Eri.fits.gz'
    x = xmrdfits(file,1)
    q=where(x.member eq 1,nmem)
    qn=where(x.member le 0)




    readcol,'data/eridanus_marla.db',  ID, ra, dec, gmag, gerr, rmag,$
            rerr, chi, sharp, rd, l, b, egr, ag, ar
    qd=where(gerr le 0.2 and rerr le 0.2 and rd le 15) 



; *******  CMD *********
     gmr=gmag-rmag
     plotsym,0,0.3,/fill
     plot,gmr[qd],rmag[qd],psym=8,yrange=[24,16],xrange=[-0.5,1.1],xstyle=1,$
          ystyle=1,/nodata,$
          xtitle='g-r',ytitle='r',title='!6 Eri: CFHT Photometry'
     loadct,0
     oplot,gmr[qd],rmag[qd],psym=8,color=80
     loadct,12

    file = '~/Dropbox/idl/mgidl/isochrones/SDSS/gircmd_t12_z0004.dat'
    readcol,file,f='f,f,f,f,f,f,f, f,f,f,f,f, f,f,f,f,f,f',$
            lg,mi,ma,ll,lr,lg,mbol,umag,gmag,rmag,imag,zmag,$
            CO,M_hec,period,pmode,logMdot,int_IMF
     b=70
     g_iso = gmag[0:b] + dmod + Ag
     r_iso = rmag[0:b] + dmod + Ar
     gr_iso=g_iso-r_iso

     c=90
     d=120
     g_isohb = gmag[c:d] + dmod + Ag
     r_isohb = rmag[c:d] + dmod + Ar
     gr_isohb=g_isohb-r_isohb


      oplot,g_iso- r_iso,r_iso
      oplot,g_isohb- r_isohb,r_isohb


     plotsym,0,0.7,/fill
     oplot,x[qn].gr,x[qn].rmag,psym=8

     plotsym,0,0.8,/fill
     oplot,x[q].gr,x[q].rmag,psym=8,color=180
     plotsym,0,0.83
     oplot,x[q].gr,x[q].rmag,psym=8

     plotsym,0,0.9,/fill
     legend,psym=[8,8],color=[180,0],['Eri Member','Non-Member'],box=0,$
            /top,/left,charsize=0.8

; *******  VHIST *********
     plothist,x[q].vcorr,bin=5,xrange=[-200,100]+100,$
     title='!6 Eri: Keck Velocities',xtitle='Velocity (km s!U-1!N)',$
     ytitle='N!Dstars!N',yrange=[0,25], ystyle=1             

     plothist,x[qn].vcorr,bin=5,/overplot
     plothist,x[q].vcorr,bin=5,/overplot,color=180,thick=2

     oplot,mw_a,mw_b*0.00014,color=100


;********** Calc Numbers for Abstract
print,'N Members = ',nmem

   mg_maxlike_sigma,x[q].vcorr,x[q].mcerr_vcorr,v,s
  ; mg_mcmc_sigma,x[q].vcorr,x[q].mcerr_vcorr,v,s,500

  r    = 16.7
  rerr = 1.1
  sig = 0.7
  serr = 0.5
  MV = -4.9
  MVerr = 0.3
  calc_wolfmass, r, rerr, sig, serr, MV,MVerr,mass,merr

stop
;******************************


if keyword_set(create_ps) then begin
   device,/close
   set_plot,'x'
endif

end
 


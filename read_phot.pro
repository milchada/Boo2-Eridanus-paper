Pro read_phot

!p.multi=[0,1,2]

readcol,'boo2_marla.db', ID, RA, DEC, g, gerr, r, rerr, chi, sharp, rd, l, b, egr, ag, ar


q=where(gerr le 0.2 and rerr le 0.2 and rd le 10) 
plot,g[q]-r[q],r[q],psym=3,yrange=[24,16],xrange=[-0.5,1.1],xstyle=1,$
          ystyle=1



readcol,'eridanus_marla.db', ID, RA, DEC, g, gerr, r, rerr, chi, sharp, rd, l, b, egr, ag, ar


q=where(gerr le 0.2 and rerr le 0.2 and rd le 10) 
plot,g[q]-r[q],r[q],psym=3,yrange=[24,16],xrange=[-0.5,1.1],xstyle=1,$
          ystyle=1
stop


end

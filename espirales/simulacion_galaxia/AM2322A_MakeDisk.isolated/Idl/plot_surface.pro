

openr,1,"surface_density.txt"
da=fltarr(5, 512)
readf,1,da
close,1


;rho0, sigma, sigmasfr, sigma_u4, sigmacheck

rho0= da(0,*)
sigma= da(1,*)
sigmacheck=da(4,*)


plot, rho0, sigma, /xlog, /ylog
oplot, rho0, sigmacheck, color=255

end


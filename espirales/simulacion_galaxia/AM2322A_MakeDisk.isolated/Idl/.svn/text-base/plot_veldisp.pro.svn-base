

openr,1,"deldisp.dat"

rsize=0L
zsize=0L
readu,1,rsize, zsize
print,"rsize=", rsize
print,"zrsize=", zsize

R=dblarr(rsize)
readu,1,R
Z=dblarr(zsize)
readu,1,Z


rho=dblarr(Zsize, Rsize)
rhophi=dblarr(Zsize, Rsize)
rhovc=dblarr(Zsize, Rsize)

readu,1,rho
readu,1,rhophi
readu,1,rhovc
kappa=dblarr(rsize)
readu,1,kappa

close,1

rho=transpose(rho)
rhophi=transpose(rhophi)
rhovc =transpose(rhovc)

rho=sqrt(rho)

plot, z, rho(0,*), xrange=[0,3]

for i=1, rsize-1 do begin

 oplot, z, rho(i,*)
  
endfor


end

 

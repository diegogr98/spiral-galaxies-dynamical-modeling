

openr,1,"gas_density.txt"

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
cumul=dblarr(Zsize, Rsize)
readu,1,rho
readu,1,cumul
close,1

rho=transpose(rho)
cumul=transpose(cumul)


openr,1,"particles.dat"
N=0L
readu,1,N
Pos=fltarr(3,N)
readu,1,Pos
close,1


dist = sqrt(pos(0,*)^2 + pos(1,*)^2)

ind=where(dist lt 0.5)

zz = pos(2,ind)


zz=zz(sort(zz))


count=indgen(n_elements(zz))+1

count= count - 0.5*max(count)

plot, z,cumul(0,*) , xrange=[0,2]

oplot, zz, count/max(double(count))*max(cumul(0,*)), linestyle=2


end

 

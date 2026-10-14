

openr,1,"gas_density.txt"
;openr,1,"gas_density_before.txt"

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
readu,1,rho
close,1

rho=transpose(rho)

plot, z, rho(1,*), xrange=[0,3]

for i=2, rsize-1 do begin

 oplot, z, rho(i,*)

endfor

;plot, r, rho(*,0), xrange=[0,13]

;oplot, r, rho(*,0), psym=5

;for i=1, zsize-1,10 do begin

; oplot, r, rho(*,i)

;endfor


end

 

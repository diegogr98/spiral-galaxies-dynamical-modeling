

openr,1,"particles.dat"
N=0L
readu,1,N
Pos=fltarr(3,N)
readu,1,Pos
close,1




plot, pos(0,0:N/2-1), pos(2,0:n/2-1), psym=3


;plot, pos(0,*), pos(2,*), psym=3




r = sqrt (pos(0,0:N/2-1)^2 +  pos(1,0:n/2-1)^2)

plot, r, pos(2,0:n/2-1) , psym=3



end

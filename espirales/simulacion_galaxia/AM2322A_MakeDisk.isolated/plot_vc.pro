
openr,1,"curve.txt"
readf,1, CC 
readf,1, V200
readf,1, LAMBDA
readf,1, MD
readf,1, JD
readf,1, MB
readf,1, DiskHeight
readf,1, BulgeSize
readf,1, R200
readf,1, H

readf,1, POINTS
R=fltarr(points)
vc=fltarr(points)
vc1=fltarr(points)
vc2=fltarr(points)
vc3=fltarr(points)
readf,1,R
readf,1,vc
readf,1,vc1
readf,1,vc2
readf,1,vc3

close,1


r=r/h

plot, r, sqrt(vc), xrange=[0,10], xstyle=1

oplot, r, sqrt(vc1), linestyle=1
oplot, r, sqrt(vc2), linestyle=2
oplot, r, sqrt(vc3), linestyle=3


end

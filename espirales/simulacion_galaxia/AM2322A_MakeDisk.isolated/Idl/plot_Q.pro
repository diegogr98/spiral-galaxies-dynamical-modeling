
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

readf,1,count
rq=fltarr(count)
Q=fltarr(count)
readf,1,rq
readf,1,Q
close,1


rq=rq/h

plot, rq, Q, xrange=[0,6],yrange=[0,4]


end

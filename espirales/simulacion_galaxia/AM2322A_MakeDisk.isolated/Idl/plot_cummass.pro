
N=4000

openr,1,"cummass.txt"
da= fltarr(5,N)
readf,1,da
close,1

r1 = da(0,*)
r2 = da(1,*)

mf = da(2,*)



plot, r1(1:*), mf(1:*), /xlog, /ylog

oplot,r2, mf, color=255

end

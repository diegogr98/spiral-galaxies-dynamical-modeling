
frun="/gmpa/mpa/vrs/Disks/disk3/"

snapbase="snap"


window,xsize=800,ysize=800


for num=0,1000,1 do begin

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname=frun+"/"+snapbase+"_"+exts
    fname=strcompress(fname,/remove_all)

;;;fname = "/gmpa/mpa/vrs/Disks/disk2/ics_disk2.dat"
    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)

    print,fname
    openr,1,fname,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
    print,"Time= ", time
    N=total(npart)
    pos=fltarr(3,N)
    readu,1,pos

    NGas=npart(0)
    NHalo=npart(1)
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)

    close,1


    if Ngas gt 0 then begin
        xgas=pos(0,0:Ngas-1)
        ygas=pos(1,0:Ngas-1)
        zgas=pos(2,0:Ngas-1)
    endif

    if Nhalo gt 0 then begin
        xhalo=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo=pos(2,0+Ngas:Nhalo+Ngas-1)
    endif

    if Ndisk gt 0 then begin
        xdisk=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    endif

    if Nbulge gt 0 then begin
        xbulge=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    endif

for cc=0,0 do begin

 r=sqrt(xhalo^2 + yhalo^2 + zhalo^2)
 ind=where(r lt 1.0)
 cx= total(xhalo(ind),/double)/n_elements(ind)
 cy= total(yhalo(ind),/double)/n_elements(ind)
 cz= total(zhalo(ind),/double)/n_elements(ind)

 print, cx, cy, cz
; xhalo = xhalo - cx
; yhalo = yhalo - cy
; zhalo = zhalo - cz
endfor
    
rmin= 0.1 
rmax= 300

BINS = 25
dlogr = alog(rmax/rmin)/BINS


 h = histogram(alog(r), min=alog(rmin), max=alog(rmax), binsize= dlogr)

r1= exp(indgen(n_Elements(h))*dlogr+alog(rmin))
r2= exp((indgen(n_Elements(h))+1.0)*dlogr+alog(rmin))
vol= 4*!PI/3.0*(r2^3-r1^3)

ind=where(h ne 0)

rho = h/vol

rho1= (h+1)/vol



x = exp((indgen(n_Elements(h))+0.5)*dlogr+alog(rmin))


   r=r(sort(r))

co= reverse(lindgen(n_elements(r)+1))


if num eq 0 then begin
  plot, x, rho, /ylog ,yrange=[min(rho1),max(rho1)],xrange=[0.1, 300],/xlog, xstyle=1
endif else begin
  oplot, x, rho, color=255
endelse


;    wait, 0.5	

endfor	

end











frun="/gmpa/mpa/vrs/Disks/disk3/"
snapbase="snap"


window,xsize=800,ysize=800


for num=0,8,1 do begin

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname=frun+"/"+snapbase+"_"+exts
    fname=strcompress(fname,/remove_all)

    
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



    xdisk=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    ydisk=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    zdisk=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)


    r=sqrt(xdisk^2+ydisk^2)

rmin= 0.1 
rmax= 40

BINS = 20
dlogr = alog(rmax/rmin)/BINS


 h = histogram(alog(r), min=alog(rmin), max=alog(rmax), binsize= dlogr)

r1= exp(indgen(n_Elements(h))*dlogr+alog(rmin))
r2= exp((indgen(n_Elements(h))+1.0)*dlogr+alog(rmin))
area= !PI*(r2^2-r1^2)

ind=where(h ne 0)

rho = h/area

rho1= (h+1)/area



x = exp((indgen(n_Elements(h))+0.5)*dlogr+alog(rmin))



if num eq 0 then begin
  plot, x, rho, /ylog ,yrange=[min(rho1),max(rho1)],xrange=[0.1, 20], xstyle=1
endif else begin
  oplot, x, rho, color=255
endelse





endfor	

end










frun="/gmpa/mpa/vrs/Disks/disk3/"
snapbase="snap"


window,xsize=800,ysize=800

xlen=17.0

for num=0,1000,1 do begin

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
   if Nbulge gt 0 then begin
        xbulge=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    endif
	  
    plot, xdisk, ydisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1

;    oplot, xbulge, zbulge, psym=3, color=255

    wait, 0.5	

endfor	

end









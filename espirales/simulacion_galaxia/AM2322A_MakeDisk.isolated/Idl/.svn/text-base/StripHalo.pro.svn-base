
fin= "../../ICs/lr_gB_eA_vc4_gf100.dat"
fout="../../ICs/lr_gB_eA_vc4_gf100_nohalo.dat"

    
npart=lonarr(6)	
  massarr=dblarr(6)
  time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartTotal=lonarr(6)	
  bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
la=intarr(bytesleft/2)

  print,fin
  openr,1,fin,/f77_unformatted
  readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
  print,"Time= ", time

  print, npart

  NGas=npart(0)
  NHalo=npart(1)
  NDisk=npart(2)
  NBulge=npart(3)
  NStars=npart(4)
  NBH=npart(5)


  pos0=fltarr(3,Ngas)
  pos1=fltarr(3,Nhalo)
  pos2=fltarr(3,Ndisk+Nbulge+nstars+NBH)

  readu,1,pos0, pos1, pos2

  vel0=fltarr(3,Ngas)
  vel1=fltarr(3,Nhalo)
  vel2=fltarr(3,Ndisk+Nbulge+nstars+NBH)

  readu,1,vel0, vel1, vel2


  id0=lonarr(Ngas)
  id1=lonarr(Nhalo)
  id2=lonarr(Ndisk+Nbulge+nstars+NBH)

  readu,1, id0, id1, id2


  ind=where((npart gt 0) and (massarr eq 0)) 
  if ind(0) ne -1 then begin
  Nwithmass= total(npart(ind))
  mass=fltarr(Nwithmass)
  endif else begin	
  Nwithmass= 0
  endelse

  if Nwithmass gt 0 then begin
  readu,1,mass
  endif

  u =fltarr(Ngas)
  readu,1,u

  close,1


  openw,1,fout,/f77_unformatted
  npartTotal(1)=0
  npart(1)=0
  massarr(1)=0


  print
  print, npart
  writeu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
  writeu,1,pos0, pos2
  writeu,1,vel0, vel2
  writeu,1, id0, id2

  if Nwithmass gt 0 then begin
  writeu,1,mass
  endif

  writeu,1,u

  close,1
 
end









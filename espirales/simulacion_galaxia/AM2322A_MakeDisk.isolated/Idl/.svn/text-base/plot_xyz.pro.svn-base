
  UnitLength_in_cm =        3.085678d21        ;  1.0 kpc 
  UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
  UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec 
  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
  UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
  UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
  UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
  GRAVITY   = 6.672d-8
  BOLTZMANN = 1.3806d-16
  PROTONMASS = 1.6726e-24
  G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
  Xh=0.76D  ; mass fraction of hydrogen
  HubbleParam= 0.7
  gamma= 5.0/3

     
    f= "test.dat"
    

    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartall=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
    la=intarr(bytesleft/2)

   
    openr,1,f,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
    print,npart,massarr
    print
    print,time,redshift
    print
    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
;    id=lonarr(N)
    readu,1  ,pos
    readu,1  ,vel
    readu,1 ; ,id
    ind=where((npart gt 0) and (massarr eq 0))
    if ind(0) ne -1 then begin
      Nm= total(npart(ind))
      mass=fltarr(Nm)	
      readu,1,mass
    endif

    NGas=npart(0)
    NHalo=npart(1)
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u
    endif
    close,1

    xgas=Pos(0,0:Ngas-1)
    ygas=Pos(1,0:Ngas-1)
    zgas=Pos(2,0:Ngas-1)

    xhalo=Pos(0,Ngas:Ngas+Nhalo-1)
    yhalo=Pos(1,Ngas:Ngas+Nhalo-1)
    zhalo=Pos(2,Ngas:Ngas+Nhalo-1)



    xdisk=Pos(0,Ngas+Nhalo-1:*)
    ydisk=Pos(1,Ngas+Nhalo-1:*)
    zdisk=Pos(2,Ngas+Nhalo-1:*)


    plot, xgas, ygas, psym=3


    r= sqrt(xgas^2 + ygas^2)


    plot, r, zgas, psym=3


    plot, xhalo, yhalo, psym=3

    plot, xdisk, ydisk, psym=3

end

















 

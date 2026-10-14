RSIZE= 51 
ZSIZE= 51 


frc=dblarr(ZSIZE, RSIZE)

openr,1, "frc_old.dat"
readu,1,frc
close,1


frc_new=dblarr(ZSIZE, RSIZE)

openr,1, "frc_new.dat"
readu,1,frc_new
close,1


plot, frc(*,3)

oplot, frc_new(*,3)


end

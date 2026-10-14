
Vvir=180.0

c=10.0


H=0.1
G=43007.1



Rvir= Vvir/(10*H)
Mvir= Vvir^3/(10*H*G)

Rs= rvir/10.0


a= rs * sqrt(2) * (alog(1+c)-c/(1+c))^(1.0/2)

x=(indgen(2000)+0.5)/2000.0 * Rvir*4



rho_nfw = mvir/(4*!PI*rs^3)  /(alog(1+c)-c/(1+c)) * 1/(x/rs) * 1/(1+x/rs)^2

rho_hq = mvir/(2*!PI) * a/x *1.0/(a+x)^3


plot, x, rho_nfw ,/xlog, /ylog

oplot, x, rho_hq , color=255


print, "a=", a
print, "rs=", rs

print,"hq-frac=", (rvir/(rvir+a))^2

end





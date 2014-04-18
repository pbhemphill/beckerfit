function beckerFit,x,p
  ; wrapper for beckerEtheo so it can be used in mpfitfun
  ; see below for which parameter is which
  lum = x
  lambda = p[0]
  w = p[1]
  esurf = p[2]
  tau = p[3]
  mass = p[4]
  radius = p[5]
  xi = p[6]
  crit = p[7]
  return,beckerEtheo(lum,lambda=lambda,esurf=esurf,w=w,tau=tau,mass=mass,radius=radius,xi=xi,crit=crit)
end

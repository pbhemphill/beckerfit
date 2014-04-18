function beckerHs,lum,mass=mass,radius=radius,xi=xi,lscale=lscale,mscale=mscale,rscale=rscale
; calculate h_s from Becker et al. (2012) - the emission height for the CRSF in the super-critical accretion regime
  if ~n_elements(mass) then mass = 1.4
  if ~n_elements(radius) then radius = 10.0
  if ~n_elements(xi) then xi=0.01
  if ~n_elements(lscale) then lscale = 1e37
  if ~n_elements(mscale) then mscale = 1.99e33
  if ~n_elements(rscale) then rscale = 1e5
  h = 2.28e3*(xi/0.01)*(mass*mscale/(1.4*1.99e33))^(-1.0)*(radius*rscale/(10.0*1e5))*(lum*lscale/1e37)
  return,h
end


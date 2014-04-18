function beckerHc,lum,lambda=lambda,mass=mass,radius=radius,esurf=esurf,tau=tau,escale=escale
; calculate h_c from Becker et al. (2012) - the emission height for the CRSF in the sub-critical accretion regime
  if ~n_elements(mass) then mass = 1.4
  if ~n_elements(radius) then radius = 10.0
  if ~n_elements(lambda) then lambda=0.1
  if ~n_elements(esurf) then esurf = 25.0
  if ~n_elements(tau) then tau = 20.0
  if ~n_elements(lscale) then lscale = 1e37
  if ~n_elements(mscale) then mscale = 1.99e33
  if ~n_elements(rscale) then rscale = 1e5
  if ~n_elements(escale) then escale = 1
  g = 6.6738480e-11 ; m^3 kg^-1 s^-2
  c = 2.99e8 ; m s^-1
  z = g*(mass*mscale/1000.D)/(c*c*(radius*rscale/100.D))
  b = (1.0+z)*esurf*escale/11.57
  h = 1.48e5*(lambda/0.1)^(-1.0)*(tau/20.0)*(mass*mscale/(1.4*1.99e33))^(19.0/14.0)*(radius*rscale/(10.0*1e5))^(1.0/14.0)*b^(-4.0/7.0)*(lum*lscale/1e37)^(-5.0/7.0)
  return,h
end


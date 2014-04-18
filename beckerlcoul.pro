function beckerLcoul,lambda=lambda,mass=mass,radius=radius,esurf=esurf,tau=tau,rscale=rscale,mscale=mscale,escale=escale
; calculate L_coul from Becker et al. (2012) - the luminosity below which
; Coulomb interactions in the plasma are alone sufficient to halt the inflow of
; plasma
  if ~n_elements(mass) then mass = 1.4
  if ~n_elements(radius) then radius = 10.0
  if ~n_elements(lambda) then lambda=0.1
  if ~n_elements(esurf) then esurf = 25.0
  if ~n_elements(tau) then tau = 20.0
  if ~n_elements(mscale) then mscale = 1.99e33
  if ~n_elements(rscale) then rscale = 1e5
  if ~n_elements(escale) then escale = 1
  g = 6.6738480e-11 ; mass^3 kg^-1 s^-2
  c = 2.99e8 ; mass s^-1
  z = g*(mass*mscale/1000.D)/(c*c*(radius*rscale/100.D))
  b = (1.0+z)*esurf*escale/11.57
  return,1.17*(lambda/0.1)^(-7.0/12.0)*(tau/20.0)^(7.0/12.0)*(mass*mscale/(1.4*1.99e33))^(11.0/8.0)*(radius*rscale/(10.0*1e5))^(-13.0/24.0)*b^(-1.0/3.0)
end


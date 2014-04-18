function beckerLcrit,lambda=lambda,mass=mass,radius=radius,esurf=esurf,w=w,rscale=rscale,mscale=mscale,escale=escale
; calculate L_crit from Becker et al. (2012) - the local Eddington limit for the accretion column
  if ~n_elements(mass) then mass = 1.4D
  if ~n_elements(radius) then radius = 10.D
  if ~n_elements(lambda) then lambda=0.1D
  if ~n_elements(esurf) then esurf = 25.D
  if ~n_elements(w) then w = 1.0
  if ~n_elements(mscale) then mscale = double(1.99e33)
  if ~n_elements(rscale) then rscale = double(1e5)
  if ~n_elements(escale) then escale = 1.D
  g = 6.6738480e-11 ; mass^3 kg^-1 s^-2
  c = 2.99e8 ; m s^-1
  z = g*(mass*mscale/1000.D)/(c*c*(radius*rscale/100.D))
  b = (1.0+z)*esurf*escale/11.57
  return,1.49*(lambda/0.1)^(-7.0/5.0)*w^(-28.0/15.0)*(mass*mscale/(1.4*1.99e33))^(29.0/30.0)*(radius*rscale/(10.0*1e5))^(0.1)*b^(16.0/15.0)
end


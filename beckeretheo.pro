function beckerEtheo,lum,lambda=lambda,esurf=esurf,w=w,tau=tau,mass=mass,radius=radius,xi=xi,crit=crit,escale=escale,mscale=mscale,rscale=rscale
;+
; NAME:
;         beckerEtheo
; PURPOSE
;         computes the theoretical cyclotron line energy for a given
;         luminosity in the Becker et al. (2012) framework, using eqns.
;         40, 51, and 58 from that work.
;
; INPUTS
;         lum: array of luminosities
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;         lambda: Lambda parameter from Becker et al. (2012). 0.1 is disk accretion, 1.0 is wind accretion.
;         w: Characterizes mean photon energy as E = wkT (w=1 is bremsstrahlung, w=3 is Planck)
;         esurf: Cyclotron line energy at surface
;         tau: Thompson optical depth required to stop infalling material
;         mass: Mass of neutron star
;         radius: Radius of neutron star
;         xi: xi parameter from Becker et al. (2012) (eqn. 36). Should be ~0.01 or therabouts.
;         escale: Scale for CRSF energy (default: 1; i.e., esurf = 30.0 is a 30 keV CRSF)
;         mscale: Scale for mass (default: 1.99e33; i.e. setting mass=1 is 1 solar mass)
;         rscale: Scale for radius (default: 1e5; i.e. r=10.0 is a 10 km radius)
;         crit: Flag to determine what gets calculated:
;               0: both super- and sub-critical accretion possible, depending on where L_crit is
;               1: only super-critical accretion
;               2: only sub-critical accretion
;
; CALLS:
;         beckerHs()
;         beckerHc()
;         beckerLcrit()
;         
; MODIFICATION HISTORY:
;         Version 0.2, 2014/04/02 Paul Hemphill (pbhemphill@physics.ucsd.edu)
;           Added w parameter to BECKERLCRIT
;         Version 0.1, 2014/03/18 Paul Hemphill (pbhemphill@physics.ucsd.edu)
;         
;
  if ~n_elements(mass) then mass = 1.4D
  if ~n_elements(radius) then radius = 10.0D
  if ~n_elements(lambda) then lambda=0.1D
  if ~n_elements(esurf) then esurf = 25.0D
  if ~n_elements(w) then w=1.0
  if ~n_elements(xi) then xi=0.01D
  if ~n_elements(tau) then tau = 20.0D
  if ~n_elements(mscale) then mscale = double(1.99e33)
  if ~n_elements(rscale) then rscale = double(1e5)
  if ~n_elements(escale) then escale = 1D
  if ~n_elements(crit) then crit = 0

  ls = beckerLcrit(lambda=lambda,mass=mass,radius=radius,esurf=esurf,w=w,rscale=rscale,mscale=mscale,escale=escale)

  if ~crit then begin ; crit=0 uses both super- and sub-critical
    sup = where(lum GE ls)
    sub = where(lum LT ls)
    e = dblarr(n_elements(lum))
    if sup[0] NE -1 and sub[0] NE -1 then begin
      hs = beckerHs(lum[sup],mass=mass,radius=radius,xi=xi,lscale=lscale,mscale=mscale,rscale=rscale)
      hc = beckerHc(lum[sub],lambda=lambda,mass=mass,radius=radius,esurf=esurf,tau=tau,escale=escale)
      e[sup] = (1.0 + hs/(radius*rscale))^(-3.0)*esurf*escale
      e[sub] = (1.0 + hc/(radius*rscale))^(-3.0)*esurf*escale
    endif else if sup[0] NE -1 then begin
      hs = beckerHs(lum,mass=mass,radius=radius,xi=xi,lscale=lscale,mscale=mscale,rscale=rscale)
      e = (1.0 + hs/(radius*rscale))^(-3.0)*esurf*escale
    endif else begin
      hc = beckerHc(lum,lambda=lambda,mass=mass,radius=radius,esurf=esurf,tau=tau,escale=escale)
      e = (1.0 + hc/(radius*rscale))^(-3.0)*esurf*escale
    endelse
  endif else if crit EQ 1 then begin ; supercritical only
      hs = beckerHs(lum,mass=mass,radius=radius,xi=xi,lscale=lscale,mscale=mscale,rscale=rscale)
      e = (1.0 + hs/(radius*rscale))^(-3.0)*esurf*escale
  endif else if crit EQ 2 then begin ; subcritical only
      hc = beckerHc(lum,lambda=lambda,mass=mass,radius=radius,esurf=esurf,tau=tau,escale=escale)
      e = (1.0 + hc/(radius*rscale))^(-3.0)*esurf*escale
  endif else begin
    print,'Error: crit must be 0, 1, or 2'
    return,'error'
  endelse
  return,e
end


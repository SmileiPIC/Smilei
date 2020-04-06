import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Physical constants in SI units
c_SI    = 299792458.
me_SI   = 9.10938356e-31
e_SI    = 1.60217662e-19
hbar_SI = 1.054571800e-34
eps0_SI = 8.854187817e-12
af      = e_SI**2/(4.*np.pi*eps0_SI*hbar_SI*c_SI)
te_SI   = e_SI**2/(4.*np.pi*eps0_SI*me_SI*c_SI**3)
Palpha  = 2./3. * af**2 / (S.namelist.wr*te_SI)
sim_hyper_volume = S.namelist.Lx

# quantum correction function g(chi) -- fit by Ridgers
def g(chi):
    return (1+4.8*(1.+chi)*np.log(1+1.7*chi)+2.44*chi**2)**(-2./3.);

# estimate the power radiated away by a single electron
g0       = S.namelist.g0
chi0     = S.namelist.chi0
n0       = S.namelist.n0
l0       = S.namelist.l0
Prad     = Palpha * chi0**2 * g(chi0)
Prad_tot = Prad * n0 * l0

# extract the power radiated away from the RadiationSpectrum diagnostics
spc_noRR_t = np.array( S.RadiationSpectrum(0).getData() )
spc_LL_t   = np.array( S.RadiationSpectrum(1).getData() )
spc_cLL_t  = np.array( S.RadiationSpectrum(2).getData() )
spc_FP_t   = np.array( S.RadiationSpectrum(3).getData() )
spc_MC_t   = np.array( S.ParticleBinning(10).getData() ) * sim_hyper_volume    ### Due to normalization of ParticleBinning (will be rediscussed w/ Fred)
gaxis      = np.array( S.ParticleBinning(10).getAxis("gamma") )

spc_noRR = np.mean(spc_noRR_t,axis=0)
spc_LL   = np.mean(spc_LL_t,  axis=0)
spc_cLL  = np.mean(spc_cLL_t, axis=0)
spc_FP   = np.mean(spc_FP_t,  axis=0)
spc_MC   = spc_MC_t[-1]/S.namelist.Tsim

dgaxis     = np.zeros(gaxis.size)
dgaxis[0]  = gaxis[0]
for i in range(1,gaxis.size): dgaxis[i] = gaxis[i]-gaxis[i-1]

integrand  = dgaxis*spc_noRR
Prad_noRR  = np.sum(integrand)
integrand  = dgaxis*spc_LL
Prad_LL    = np.sum(integrand)
integrand  = dgaxis*spc_cLL
Prad_cLL   = np.sum(integrand)
integrand  = dgaxis*spc_FP
Prad_FP    = np.sum(integrand)

# extract the power radiated away from the MC ParticleBinning diagnostics

integrand  = dgaxis*spc_MC
Prad_MC    = np.sum(integrand)

# Physical constants in SI units
c_SI    = 299792458.
me_SI   = 9.10938356e-31
e_SI    = 1.60217662e-19
hbar_SI = 1.054571800e-34
eps0_SI = 8.854187817e-12
af      = e_SI**2/(4.*np.pi*eps0_SI*hbar_SI*c_SI)
te_SI   = e_SI**2/(4.*np.pi*eps0_SI*me_SI*c_SI**3)
Palpha  = 2./3. * af**2 / (S.namelist.wr*te_SI)

# quantum correction function g(chi) -- fit by Ridgers
def g(chi):
    return (1+4.8*(1.+chi)*np.log(1+1.7*chi)+2.44*chi**2)**(-2./3.);

# estimate the power radiated away by a single electron
g0       = S.namelist.g0
chi0     = S.namelist.chi0
n0       = S.namelist.n0
l0       = S.namelist.l0
Prad     = Palpha * chi0**2 * g(chi0)
Prad_tot = Prad * n0 * l0

# extract the power radiated away from the RadiationSpectrum diagnostics
spc_noRR_t = np.array( S.RadiationSpectrum(0).getData() )
spc_LL_t   = np.array( S.RadiationSpectrum(1).getData() )
spc_cLL_t  = np.array( S.RadiationSpectrum(2).getData() )
spc_FP_t   = np.array( S.RadiationSpectrum(3).getData() )
spc_MC_t   = np.array( S.ParticleBinning(10).getData() ) * sim_hyper_volume    ### Due to normalization of ParticleBinning (will be rediscussed w/ Fred)
gaxis      = np.array( S.ParticleBinning(10).getAxis("gamma") )

spc_noRR = np.mean(spc_noRR_t,axis=0)
spc_LL   = np.mean(spc_LL_t,  axis=0)
spc_cLL  = np.mean(spc_cLL_t, axis=0)
spc_FP   = np.mean(spc_FP_t,  axis=0)
spc_MC   = spc_MC_t[-1]/S.namelist.Tsim

dgaxis     = np.zeros(gaxis.size)
dgaxis[0]  = gaxis[0]
for i in range(1,gaxis.size): dgaxis[i] = gaxis[i]-gaxis[i-1]

integrand  = dgaxis*spc_noRR
Prad_noRR  = np.sum(integrand)
integrand  = dgaxis*spc_LL
Prad_LL    = np.sum(integrand)
integrand  = dgaxis*spc_cLL
Prad_cLL   = np.sum(integrand)
integrand  = dgaxis*spc_FP
Prad_FP    = np.sum(integrand)

# extract the power radiated away from the MC ParticleBinning diagnostics

integrand  = dgaxis*spc_MC
Prad_MC    = np.sum(integrand)

# Compare
print( "chi0          = ", chi0)
print( "g(chi0)       = ", g(chi0))
print( "Prad (theory) = ", Prad_tot)
print( "Prad (noRR)   = ", Prad_noRR)
print( "Prad (LL)     = ", Prad_LL)
print( "Prad (cLL)    = ", Prad_cLL)
print( "Prad (FP)     = ", Prad_FP)
print( "Prad (MC)     = ", Prad_MC)

Validate("Prad (noRR)",Prad_noRR, 1.e-6*Prad_noRR)
Validate("Prad   (LL)",  Prad_LL, 1.e-6*Prad_LL  )
Validate("Prad  (cLL)", Prad_cLL, 1.e-6*Prad_cLL )
Validate("Prad   (FP)",  Prad_FP, 1.e-2*Prad_FP  )
Validate("Prad   (MC)",  Prad_MC, 1.e-2*Prad_MC  )

# Table to check
Validate("Radiation Spectrum (noRR)", spc_noRR, 1.e-6)
#spc_LL
#spc_cLL
#spc_FP
#spc_MC

import happi
from numpy import array, log10

S = happi.Open(["./restart*"], verbose=False)

v = []
N = []
Nref = []
Npart = []
for i in range(0,len(S.namelist.DiagParticleBinning),2):
    v    += [S.namelist.Species["Da_%d"%(i//2)].mean_velocity[0]]
    N    += [S.ParticleBinning(i  ,sum={"x":"all"}).getData()[-1]]
    Nref += [S.ParticleBinning(i+1,sum={"x":"all"}).getData()[-1]]
    Npart += [S.Scalar("Ntot_He_%d"%(i//2)).getData()[-1]]

Validate("log10 of created helium3 density", log10(array(N)+1e-40), 0.2)
Validate("number of created helium3", Npart, 50.)


# # Some ploting to compare w theory
# E_theory = array([
#     1.00000000e+03, 1.49328996e+03, 2.22991489e+03, 3.32990952e+03,
#     4.97252044e+03, 7.42541482e+03, 1.10882974e+04, 1.65580431e+04,
#     2.47259595e+04, 3.69230270e+04, 5.51367853e+04, 8.23352077e+04,
#     1.22950339e+05, 1.83600506e+05, 2.74168792e+05, 4.09413503e+05,
#     6.11373072e+05, 9.12957268e+05, 1.36330992e+06, 2.03581701e+06,
#     3.04006509e+06, 4.53969867e+06, 6.77908643e+06, 1.01231417e+07,
#     1.51167858e+07, 2.25737444e+07, 3.37091458e+07, 5.03375289e+07,
#     7.51685263e+07, 1.12248405e+08, 1.67619416e+08, 2.50304391e+08,
#     3.73777033e+08, 5.58157490e+08
# ]) * 2.0141 # in eV
# crosssection_theory = array([
#     1.37959339e-09, 2.71808933e-07, 1.93984386e-05, 5.96656956e-04,
#     9.09434767e-03, 7.80510801e-02, 4.23198290e-01, 1.59539840e+00,
#     4.50717898e+00, 1.01191544e+01, 1.89534351e+01, 3.08794670e+01,
#     4.52831325e+01, 6.12170894e+01, 7.72386840e+01, 9.12459059e+01,
#     1.00937123e+02, 1.04962615e+02, 1.03790887e+02, 9.92555901e+01,
#     9.30691435e+01, 8.55168643e+01, 7.52550684e+01, 6.07107292e+01,
#     4.26722519e+01, 2.52075010e+01, 1.24955923e+01, 5.39117181e+00,
#     2.13561911e+00, 8.02190038e-01, 2.81582182e-01, 8.91232341e-02,
#     2.57170089e-02, 6.56861761e-03
# ]) * 1e-31 # in m^2
# mass_deuterium = 3671.46*511000.
# Wr = S.namelist.Main.reference_angular_frequency_SI # in Hz
# density =  S.namelist.Species[0].number_density * 0.0003142 * Wr**2 # in m^-3
# time = S.namelist.Main.simulation_time / Wr # in s
# vrel = 3e8 * sqrt(1. - (E_theory/mass_deuterium + 1.)**-2) # in m/s
# N_theory = vrel * density * crosssection_theory * time
# E = (1./sqrt(1.-array(v)**2)-1.) * mass_deuterium
# clf()
# loglog(E_theory,N_theory,'-')
# loglog(E, array(N) / array(Nref), '+')

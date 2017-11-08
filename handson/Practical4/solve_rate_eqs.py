##################################################################
### SOLVING SIMPLE RATE EQUATIONS FOR THE IONIZATION OF CARBON ###
##################################################################

# control parameter
tiny = 1.e-18

# conversion factor (between code units & atomic units)
au_to_w0 = 4.134137172e+16 / reference_angular_frequency_SI;
Ec_to_au = 3.314742578e-15 * reference_angular_frequency_SI;
Eau      = aL*Ec_to_au

# PIC time-teps
dtPIC = Lmu * 1.e-6/3.e8 * 1.0e15 / rest
print '********** '
print '- [theory] dtPIC = ',dtPIC,' fs'
print ' '

# Carbon atom properties
Ip    = np.zeros(Zat)
l     = np.zeros(Zat)

nstar = np.zeros(Zat)
alpha = np.zeros(Zat)
beta  = np.zeros(Zat)
gamma = np.zeros(Zat)
Wadk  = np.zeros(Zat)

Ip  = [11.2602/27.2114,24.3845/27.2114,47.8877/27.2114,64.4935/27.2114,392.0905/27.2114,489.9931/27.2114]
l   = [1,1,0,0,0,0]

for Z in range(0,Zat):
    nstar[Z] = (Z+1.0) * m.sqrt(1./2./Ip[Z])
    cst      = 2.*nstar[Z]
    alpha[Z] = cst - 1.
    beta[Z]  = 2.**(cst-1.)/cst/m.gamma(cst) * (8.*l[Z]+4.) * Ip[Z] * au_to_w0
    gamma[Z] = 2.*(2.*Ip[Z])**1.5
    Wadk[Z]  = m.sqrt(6./m.pi) * beta[Z] * (gamma[Z]/Eau)**(cst-1.5) * m.exp(-1./3. * gamma[Z]/Eau)

# Preparing arrays
t    = np.zeros(nt)
E    = np.zeros(nt)
n    = np.zeros([Zat+1,nt]); n[0,0]=1.
Wint = np.zeros(Zat)
Env  = np.zeros(nt)

# Solving the rate equations numerically
for it in range(1,nt):
    t[it]   = it*dt
    E[it]   = aL*m.sin(t[it]) * laser_time_envelope(t[it])
    Env[it] = laser_time_envelope(t[it])

    # neutral atom
    delta  = gamma[0]/( abs(E[it])*Ec_to_au)
    if (delta>tiny):
        W        = beta[0] * (delta)**alpha[0] * m.exp(-delta/3.)
        Wint[0] += W
        n[0,it]  = (1.-W*dt/2.)/(1.+W*dt/2.) * n[0,it-1]

    # from charge 1 to charge Zat-1
    if (Zat!=1):
        for Z in range(1,Zat):
            deltam    = gamma[Z-1]/( abs(E[it])*Ec_to_au) 
            deltap    = gamma[Z]  /( abs(E[it])*Ec_to_au) 
            if (deltam>tiny) and (deltap>tiny):
                Wm       = beta[Z-1] * (deltam)**alpha[Z-1] * m.exp(-deltam/3.)
                Wp       = beta[Z  ] * (deltap)**alpha[Z  ] * m.exp(-deltap/3.)
                Wint[Z] += Wp
                n[Z,it]  = (1.-Wp*dt/2.)/(1.+Wp*dt/2.)*n[Z,it-1] + Wm*dt/(1.+Wm*dt/2.)*n[Z-1,it-1]

    # last charge state
    delta = gamma[Zat-1]/( abs(E[it])*Ec_to_au)
    if (delta>tiny):
        W          = beta[Zat-1] * (delta)**alpha[Zat-1] * m.exp(-delta/3.)
        n[Zat,it]  = n[Zat,it-1]+W*dt/(1.+W*dt/2.)*n[Zat-1,it-1]

# Compare ionisation rates
Wint = Wint/nt
for Z in range(Zat):
    print '- [theory] Z =',Z,'->',Z+1,'     Wadk=',Wadk[Z]* reference_angular_frequency_SI,'    W=',Wint[Z] * reference_angular_frequency_SI

nsum = 0.
for Z in range(Zat+1):
    nsum += n[Z,nt-1]
print ' '
print '- [theory] test cons. part. nb:', nsum
print '********** '



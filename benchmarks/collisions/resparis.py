import numpy as np
from math import pi as pi
import matplotlib.pyplot as plt

def resparis(te,ti,znuc,atwt,rho,p1,p2,p3,pot,g0,g1,alpha,visu,sty,numfig):
	# Ce M-file calcule la resistivite electrique utilisee dans le code hybride Paris 
	# Ref. Lee & More, POF (1984)
	# Entrees :
	# . te(1:n),ti(1:n) : temperatures electronique et ionique (keV)
	# . znuc        : numero atomique
	# . atwt        : masse atomique
	# . rho         : densite (g/cm^3)
	# . p1,p2,p3    : parametres d'ajustement
	# . pot         : potentiel d'ionisation (eV)
	# . g0,g1       : poids statistique des atomes neutre et ionise
	#                 (=(2,1) pour Al)
	# . alpha       : facteur intervenant dans l'ionisation par pression
	#                 (=1.5 pour Al)
	# . visu        : affichage des resultats si > 0
	# . sty         : style du trace ('-','--',':','-.'...) si visu >0
	# Sorties :
	# . eta         : resistivite (Ohm.m)
	# . teff        : temps moyen de collision e-i (s)
	# . ze          : degre moyen d'ionisation
	# . xne         : densite electronique (1/cm^3)
	# . efev        : energie de Fermi (eV)
	#
	# Caracteristiques de Al
	# atwt=27;
	# znuc=13;
	# rho=2.7;
	# p1=1;
	# p2=1;
	# p3=0.46;
	# pot =5.98 eV
	# (g0,g1) =(2,1)
	# alpha =1.5
	# Caracteristiques de Cu
	# atwt=63.5;
	# znuc=29;
	# rho=8.96;
	# p1=1;
	# p2=1;
	# p3=0.11;
	# pot = 7.73d0
	# Caracteristiques de Au
	# atwt=197;
	# znuc=79;
	# rho=19.3;
	# p1=1;
	# p2=1;
	# p3=0.07;
	# temin,temax : bornes min et max de la temperature (keV)

	q = 4.8032e-10 # e/sqrt(4 pi e0) in cm^1.5*g^0.5/s
	m = 9.1094e-28

	eta = 0*te; teff = 0*te; ze = 0*te; xne = 0*te; efev = 0*te

	n = len(te)
	if len(ti) !=n:
		print('Erreur : len(te) != len(ti) !')
		return

	ze = np.zeros(n)
	xne = np.zeros(n)
	zlam = np.zeros(n)
	teff = np.zeros(n)
	efev = np.zeros(n)
	aa = np.zeros(n)
	eta = np.zeros(n)
	wp = np.zeros(n)
	if visu > 0:
		print ("{:^10s} "*12).format('Te','ne','fe','ztf','ze','teff','efev','eta','iflag','zeta','zlam','aa')
	for i in range(n):
		# Calcul de l'ionisation
		ztf = TFMore(znuc,atwt,rho,te[i])
		fe,a = Saha(znuc,atwt,rho,te[i],pot,g0,g1,alpha)
		ze[i] = fe**(2./ztf**2)*ztf + (1.-fe**(2./ztf**2))*fe
		# ze[i] = ztf
		# ze[i] = 1.
		# Resistivite en ohm.m     
		zeta,xne[i],zlam[i],teff[i],efev[i],iflag = relaxt(te[i],ti[i],znuc,atwt,rho,ze[i],p1,p2,p3)
		cond,aa[i] = sigc(ze[i],teff[i],zeta,xne[i])
		eta[i] = 1./cond
		wp[i] = (4.*pi*xne[i]*q**2/m)**0.5
		if visu > 0:
			print ("{:10.4g} "*12).format(te[i],xne[i],fe,ztf,ze[i],teff[i],efev[i],eta[i],iflag,zeta,zlam[i],aa[i])

	if visu:
	
		f = plt.figure(numfig)
		ax = f.add_subplot(1,1,1)
		ax.loglog(te,eta,sty)
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$\eta\;(\Omega\cdot m)$")
		f.canvas.draw()

		f = plt.figure(numfig+1)
		ax = f.add_subplot(1,1,1)
		ax.loglog(te,1./eta,sty)
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$\sigma\;(\Omega^{-1}\cdot m^{-1})$")
		f.canvas.draw()

		f = plt.figure(numfig+2)
		ax = f.add_subplot(1,1,1)
		ax.loglog(te,1./teff,sty)
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$\nu_{ei}\;(s^{-1})$")
		f.canvas.draw()

		f = plt.figure(numfig+3)
		ax = f.add_subplot(1,1,1)
		ax.loglog(te,wp*teff,sty)
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$\omega_{p}\tau_{ei}$")
		f.canvas.draw()

		f = plt.figure(numfig+4)
		ax = f.add_subplot(1,1,1)
		ax.loglog(te,np.maximum(teff,8.85e-12*eta),sty)
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$\tau_{n}$ (s)")
		f.canvas.draw()

		f = plt.figure(numfig+5)
		ax = f.add_subplot(1,1,1)
		ax.plot(te,ze,sty);
		ax.set_xlim(te.min(),te.max())
		ax.set_xlabel('T (keV)')
		ax.set_ylabel(r"$Z^{*}$")
		f.canvas.draw()
		
		plt.show()
	
	return eta,teff,ze,xne,efev,aa,zlam


def sigc(ze,teff,zeta,xne):
	#...................................................................
	#  this subroutine calculates electrical conductivity 1/(ohm.m).
	#  the program was written by yim t. lee h-div. llnl april ,1983 .
	#  revised by rmm, 0   1/87
	#   output:
	#    sigc  = electrical conductivity, 1/(ohm*m)
	#    sigc is a function of:
	#         teff = relaxation time (sec)
	#          xne = electron density (/cc)
	#         zeta = degeneracy parameter = ln(1+exp(mu/kT))
	#           ze = ionization state
	#                 teff, xne, zeta are calculated in subr. relaxt
	#                 ze is calculated in subr. zcalc
	#  coefficient [e**2/m = 2.8179e-4 per (ohm*cm)]
	#  aa should vary between 1 (degenerate) and 32/3pi
	#  verified to < 1% error of the actual integral; rmm 12/86
	#                                                 see tknchk
	#...................................................................
	aa   = ((0.129*zeta+0.347)*zeta+3.39)/((0.124*zeta+0.511)*zeta+1.0)
	sig  = 2.8179e-4*xne*teff*aa
	rzb = 1./(ze*(1.0+zeta))
	# Prise en compte heuristique des collisions e-e (preferer la correction de Decoster dans relaxt)
	z = 100.*sig/(1.0+rzb) 
	return z,aa

	
def Saha(znuc,atwt,rho,tek,pot,g0,g1,alpha):
	# Ionisation deduite du modele de Saha
	# Ref. : M. Desjarlais, Contrib. Plasma Phys. 41, 267 (2001)
	# Entrees :
	# . te    : temperature (keV)
	# . rho   : densite (g/cm^3)
	# . atwt  : masse atomique
	# . znuc  : numero atomique
	
	te = tek*1e3
	na = (rho/atwt)*6.0221e23
	ra = (3./(4.*pi*na))**(1./3.)
	a = (alpha*1.4399e-7/(pot*ra))**1.5
	b = np.exp(-(pot/te)*(1.-a))
	k = 6.0372e+21*(g1/g0)*(te**1.5/na)*b
	if k<1e10:
		fe = 0.5*((k**2+4.*k)**0.5-k)
	else:
		fe = 1.-1./k
	return fe,a

def relaxt(temk,tik,znuc,atwt,rho,ze,p1,p2,p3):
	#...................................................................
	#  this subroutine calculates electron relaxation time
	#  equation references to Lee & More, Phys. Fluids 27, 1273 (1984).
	#
	#  inputs:
	#    znuc = nuclear charge
	#    atwt = atomic weight
	#    temk = temperature in kev
	#    rho  = mass density (gram/cm**3)
	#    ze   = ionization state (calculated in subr. zcalc)
	#    tik   = ion temperature, kev;  ( tik  = tek for now;  01/13/87)
	#  output:
	#          zeta   = degeneracy parameter = ln(1+exp(mu/kT))
	#          xne    = electron density/cc
	#          zlam   = Coulomb logarithm
	#          rlxtem = electron relaxation time (sec)
	#          iflag specifies the region:  1-5:
	#...................................................................
	#  esqr is in (ev-A):
	#  p1 is a multiplier for the minimum screening length:
	#  p2 controls the minimum mean free path:
	#  p3 modifies the Bloch-Gruneisen formula:
	#  iflag specifies the region:  1-5
	#  iflag = 0 is neutral dominated (ze < 1)
	#  iflag = 1 is Spitzer plasma
	#  iflag = 2 for strongly coupled (bmax=ri)
	#  iflag = 3 when log lambda = 2
	#  iflag = 4 at mfp=ri
	#  iflag = 5 for liquid metal
	#  iflag = 6 for solid metal
	amu =1.66042
	esqr = 14.39964 # e^2/(4*pi*e0) in eV*angstrom
	angs=1.e-8
	iflag = 1
	# temperatures (eV)
	te = 1.e3*temk
	ti = 1.e3*tik
	#ti=1000/1.1604e4
	
	# densite ionique      
	xni = (rho/(amu*atwt))*1e24
	# densite electronique
	xne = ze*xni
	# distance interatomique      
	ri    = (3./(4.*pi*xni))**(1./3.)
	# energie de Fermi (eV)    
	efev  = 3.81004*(3.*pi*pi*xne*(1.e-24))**(2./3.)
	x     = (efev/te)**0.5
	
	#   zimmerman analytic fit to ln(1.+exp(mu/kT)):
	zeta=(0.7531+x*(0.1679+x*0.3108))*x**3/(1.0+x*(0.2696+x*(0.2280+x*0.3099)))
	#   form screening length as in Eq. (19):
	#   2/3 in ef term corrects degenerate (Mott) limit.
	#   aee = 1/(electron screening length(cm))**2
	aee = 4.*pi*esqr*angs*xne/np.sqrt(te**2+(2.*efev/3.)**2)
	#   aii = 1/(ion debye length(cm))**2
	#   formula for aii gives correct ion Debye length even when (ze < 1)
	ze1 = max(ze,1.)
	aii = 4.*pi*ze1*esqr*angs*xne/ti
	#   form "total" screening length (cm)
	bmax=1.0/np.sqrt(aee+aii)
	#   minimum screening length is the ion-sphere radius*p1
	if bmax < ri*p1: iflag = 2
	bmax   = max(bmax,ri*p1)
	#   quantum and classical formulas for minimum impact parameter
	#   see Eq. (20-22)
	bminq  = angs*(3.81004*(4.*pi*pi/3.)/te)**0.5
	#bmincl = ze1*angs*esqr/(3.*te)
	bmincl = ze*angs*esqr/(3.*te)
	bmin   = (bminq**2 + bmincl**2)**0.5
	#   zlam is the logarithm of lamda  Eq. (17):
	zlam = 0.5*np.log(1.0+(bmax/bmin)**2)
	#zlam=5.
	#   Coulomb log is required to be >= 2
	if zlam <= 2.: iflag = 3
	zlam = max(zlam,2.)
	#   analytic fit to (1.+exp(-mu/kT))*F(1/2)(-mu/kT) as used in Eq. (24):
	fam = (0.882-0.160*zeta)+(0.200+0.671*zeta)*zeta**0.5
	#   a0 = Bohr radius(A), ve=sqrt(kT/m) in cm/sec, e0=esqr/a0 (ev)
	a0 = 0.529177
	e0 = esqr/a0
	ve = 2.1877e+8*(te/e0)**0.5
	c1 = 3./(2.*pi*2.**0.5)
	if ze < 0.95: 
	    # to calculate electron-neutral scattering
	    # tabu1=electron-ion tau
	    # tabu2=electron-neutral tau
	    # sigma=electron-neutral cross section
	    iflag = 0
	    #   xi = ion density, xn0 = neutral density
	    xi   = xne
	    xn0  = (0.95-ze)*xni
	    xn0  = max(xn0,1.e-10)
	    vth  = ve*np.sqrt(3.0)
	    sigma =2.0e-15
	    tau2 = 1./(xn0*vth*sigma)
	    c2   = fam/zlam
	    tau1 = c1*((te/e0)**2)*(a0*angs/ve)/(xi*((angs*a0)**3))*c2
	    tabu =tau1/(1.0+tau1/tau2)
	else:
	    #   Eq. (24) for the relaxation time:
	    c2 = fam/(zlam*(ze**2))
	    tabu = c1*((te/e0)**2)*(a0*angs/ve)/(xni*((angs*a0)**3))*c2
	
	#ebar = np.sqrt((1.5*te)**2 + efev**2)
	ebar = 1.5*te
	vbar=np.sqrt(ebar/13.605)*(2.1877e8)
	# minimum mean free path is p2*ion-sphere radius:
	tmin= p2*ri/vbar
	tb  =max(tmin,tabu)
	if tmin > tabu: iflag = 4
	
	# Contribution e-e dans le regime plasma
	# alpha0 : parametre correctif (A. Decoster, Modelling of collisions)
	#ztest=1
	#zion=max(ze,ztest)
	#alpha0=((zion**2)+2.9821*zion+1.0822)/(zion**2+1.603*zion+0.904)      
	#taueep=tb/(alpha0-1.)
	# Temps de collision plasma total
	#tb = tb*taueep/(tb+taueep)
	
	# # Cowan melting temperature, tmelt (ev):
	# bc=0.6*znuc**(1.0/9.0)
	# xsi=9.0*(znuc**0.3)*(amu*1.e-24)*xne/ze
	# xsf=(xsi/(1.0+xsi))**4
	# tmelt=xsf*0.32*xsi**(2.0*bc-2/3.)
	# # Bloch Gruneisen law: Eq. (32)
	# # p3 helps tune this to match known solid resistivity
	# tmetal=p3*50.0*(ri/vbar)*tmelt/ti
	# if tmetal > tb:
	# 	iflag = 6
	# 	if ti > tmelt:
	# 		iflag = 5
	# # above melting temperature:  (Eq. 33)
	# if ti >tmelt:
	# 	tmetal=tmetal/1.35
	# 
	# # Calcul de tau-ee metal (formule de Fisher)
	# hbar=1.0546e-34    # hbar en J.s
	# akb=1.3807e-23		# cste de Boltzmann en J/K
	# tek=1.1604d4*te    # tek en K
	# efevj=efev*1.6022e-19
	# aee2=((akb*tek)**2)/(hbar*efevj)
	# taueem=1./aee2
	# # Temps de collision metal total
	# tmetal = tmetal*taueem/(tmetal+taueem)
	# 
	# rlxtem = max(tmetal,tb)
	
	rlxtem = tb
	
	return zeta,xne,zlam,rlxtem,efev,iflag



def TFMore(Z,A,rho,Te):
	# Cette routine calcule l'etat d'ionisation Zeff a partir du
	# fit de Thomas-Fermi donne par More.
	# Te (keV)
	# rho (g/cm^3)
	
	alpha=14.3139
	beta=0.6624
	a1=0.003323
	a2=0.971832
	a3=9.26148e-5
	a4=3.10165
	b0=-1.7630
	b1=1.43175
	b2=0.315463
	c1=-0.366667
	c2=0.983333
	
	T0=1000.*Te/Z**(4./3.)
	R=rho/(Z*A)
	Tf=T0/(1.+T0)
	A=a1*(T0**a2)+a3*(T0**a4)
	B=-np.exp(b0+b1*Tf+b2*Tf**7)
	C=c1*Tf+c2
	Q1=A*(R**B)
	Q=(R**C+Q1**C)**(1./C)
	x=alpha*(Q**beta)
	Zeff=Z*x/(1.+x+np.sqrt(1.+2.*x))
	
	return Zeff


execfile("scripts/Diagnostics.py")

S = Smilei("tests/tst6_diagnostics")

A = S.ParticleDiagnostic(0, timesteps=[20,80])
A.plot(data_min=0, data_max=1000)

S.ParticleDiagnostic(0, timesteps=0, slice={"x":[3,4]} ).plot(figure=1)
S.ParticleDiagnostic(0, timesteps=0, slice={"vx":"all"}).plot(figure=2, xmax=1, vmin=0, vmax=11)

A = S.ParticleDiagnostic( 0 )
B = S.ParticleDiagnostic( 1 )
multiPlot(A, B, shape=[1,2])

A = S.ParticleDiagnostic(0,slice={"x":"all"})
B = S.ParticleDiagnostic(1,slice={"x":"all"})
multiPlot(A, B)

A=S.ParticleDiagnostic(3, slice={"ekin":"all"})
B=S.ParticleDiagnostic(3, slice={"ekin":[0,0.001]})
multiPlot(A,B)

S.ParticleDiagnostic("#2/#0").plot(figure=1)
S.ParticleDiagnostic("#2/#0", slice={"x":"all","vx":"all"}).plot(figure=2)

A = S.Scalar("Etot")
A.plot()

F=S.Field("Ex")
F.plot()

S.Probe(0,"Ex").plot()


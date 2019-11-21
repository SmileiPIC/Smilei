import happi

S = happi.Open(["./restart*"], verbose=False)

Validate("ion Pyy vs x", S.ParticleBinning(0, sum={"y":"all"}).getData()[-1], 0.4e-10)
Validate("ion Pyy vs y", S.ParticleBinning(0, sum={"x":"all"}).getData()[-1], 0.4e-10)


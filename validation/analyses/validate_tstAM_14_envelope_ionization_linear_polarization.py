import happi

S = happi.Open(["./restart*"], verbose=False)

Validate("electron px", S.ParticleBinning(0, timesteps=2000).getData()[-1], 5e-10)
Validate("electron py", S.ParticleBinning(1, timesteps=2000).getData()[-1], 5e-10)
Validate("electron pz", S.ParticleBinning(2, timesteps=2000).getData()[-1], 5e-10)

import happi

def Avg(an_array):
    return sum(an_array) / len(an_array)

S = happi.Open(["./restart*"], verbose=False)

Ukin = S.Scalar("Ukin").getData()
Ukin_electron = S.Scalar("Ukin_electron").getData()
Ukin_ion =  S.Scalar("Ukin_ion").getData()
Uelm = S.Scalar("Uelm").getData()
PoyXmin = S.Scalar("PoyXmin").getData()
Utot = S.Scalar("Utot").getData()

Validate("Ukin: ", Ukin / Avg(Ukin), 0.003)
Validate("Ukin_electron: ", Ukin_electron / Avg(Ukin_electron), 0.005)
Validate("Ukin_ion: ", Ukin_ion / Avg(Ukin_ion), 0.003)
Validate("Uelm: ", Uelm / Avg(Uelm), 0.015) # varies a lot
Validate("PoyXmin: ", PoyXmin / Avg(PoyXmin), 0.0015) # varies a lot
Validate("Utot: ", Utot / Avg(Utot), 0.005)

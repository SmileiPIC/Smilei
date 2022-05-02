import happi

def Avg(an_array):
    return sum(an_array) / len(an_array)

S = happi.Open(["./restart*"], verbose=False)

# from pprint import pprint
# pprint(vars(S.namelist.DiagScreen))

# Scalars, M. Lobet recommended to first check the energy
ukin = S.Scalar("Ukin").getData()
uelm = S.Scalar("Uelm").getData()
utot = S.Scalar("Utot").getData()

Validate("Ukinetic energy evolution: ", ukin / Avg(ukin), 1e-3)
Validate("Uelectromag evolution: ", uelm / Avg(uelm), 0.006)
Validate("Total energy evolution: ", utot / Avg(utot), 1e-3)

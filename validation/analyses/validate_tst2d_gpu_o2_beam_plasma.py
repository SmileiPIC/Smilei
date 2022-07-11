import happi


def Avg(an_array):
    return sum(an_array) / len(an_array)


S = happi.Open(["./restart*"], verbose=False)

a_scalar_list = ["Utot", "Uelm", "Ukin", "Uelm_Ex", "Uelm_Ey", "Uelm_Ez", "Uelm_Bx_m",
                 "Uelm_By_m", "Uelm_Bz_m", "Ukin_electron-beam", "Ukin_electron", "Ukin_ion"]
a_scalar_tolerance = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
                      0.003, 0.003, 0.003, 0.003, 0.003]

for i, a_scalar in enumerate(a_scalar_list):
    scalar_data = S.Scalar(a_scalar).getData()
    Validate(a_scalar, scalar_data / Avg(scalar_data), a_scalar_tolerance[i])

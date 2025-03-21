from RecombinationHistory import RecombinationHistory

Yp = 0.245
z_reion = 8
delta_z_reion = 0.5
z_Hereion = 3.5
delta_z_Hereion = 0.5

rec = RecombinationHistory(Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion)
rec.solve()

# find when tau = 1
x_tau = 0
tau = 0
while tau < 1:
    x_tau += 0.01
    tau = rec.tau_of_x(x_tau)

print(f"x_tau = {x_tau}")

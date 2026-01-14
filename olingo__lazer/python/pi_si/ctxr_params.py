from math import sqrt, exp
from sage.all import var, solve

vname = "param"

secpar = 256

if secpar == 128:
    ringdim = 256
    q = 81623040078337
    q_hat = 23473657271235664258052971185590273
    moddim = 8
    secdim = 8
    sigma_y = 2^10
    sigma_w = 2^29
    kappa_y = 29
    kappa_w = 32
    nu = 23
    moddim_hat = secdim_hat = m = 21
    sigma_tdec = 123298420299561636
elif secpar == 192:
    ringdim = 256
    q = 174243429534721
    q_hat = 115233010602376640884019849787898369
    moddim = 11
    secdim = 11
    sigma_y = 2^10
    sigma_w = 2^30
    kappa_y = 29
    kappa_w = 33
    nu = 39
    moddim_hat = secdim_hat = m = 29
    sigma_tdec = 283537181095073949
elif secpar == 256:
    ringdim = 256
    q = 254716596705793
    q_hat = 327595548451764648065623801250983937
    moddim = 14
    secdim = 14
    sigma_y = 2^10
    sigma_w = 2^30
    kappa_y = 29
    kappa_w = 33
    nu = 60
    moddim_hat = secdim_hat = m = 37
    sigma_tdec = 551404451102749325

rej_rate = 0.01
tail_bound_2 = 1.01
while ((tail_bound_2**ringdim) * exp((ringdim / 2) * (1 - tail_bound_2**2))) > rej_rate:
    tail_bound_2 += 0.01


deg   = 256  # ring Rp degree d
mod   = q_hat    # ring Rp modulus p
dim   = (secdim_hat + m, moddim_hat + secdim_hat + m + 2*secdim)  # dimensions of A in Rp^(m,n)

# Partition of s
wpart = (
    [[i] for i in range(moddim_hat + secdim_hat + m)] +
    [[i] for i in range(moddim_hat + secdim_hat + m, moddim_hat + secdim_hat + m + secdim)] + # first 6 elements of m0
    [[i] for i in range(moddim_hat + secdim_hat + m + secdim, moddim_hat + secdim_hat + m + 2*secdim)] # first 6 elements of m1
)

B = tail_bound_2 * sigma_y * sqrt(ringdim)
x = var('x')
eq = x**2 *sqrt(ringdim)+sqrt(ringdim)*x - B == 0
solutions = solve(eq, x)
beta = int(solutions[0].rhs().n() if solutions[0].rhs().n() > 0 else solutions[1].rhs().n())

# l2-norm bounds
wl2 = (
    [tail_bound_2 * sqrt(2/3) * sqrt(ringdim)*4 for _ in range(moddim_hat + secdim_hat + m)] +
    [beta * sqrt(ringdim)*4 for _ in range(moddim_hat + secdim_hat + m, moddim_hat + secdim_hat + m + secdim)] +
    [beta * sqrt(ringdim)*4 for _ in range(moddim_hat + secdim_hat + m + secdim, moddim_hat + secdim_hat + m + 2*secdim)]
)


# binary coefficients (0 = no restriction)
wbin = [0 for _ in range(moddim_hat + secdim_hat + m + 2*secdim)]

# rejection sampling flags (0 = not used)
wrej = [0 for _ in range(moddim_hat + secdim_hat + m + 2*secdim)]

# Optional: linf-norm bound on s (commented out)
# wlinf = 1


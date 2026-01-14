import os
#os.chdir("/Users/phmbp2025/Desktop/lwe-estimator")

from estimator import *
from estimator.lwe_dual import dual_hybrid

# Security parameter

secpar = 256

print("Security level: ", secpar)

# Choose priority for (1) small ver key (2) small sig or (3) 'balanced' i.e. small |pk+sig|

balanced = 1
small_pk = 0
small_sig = 0 


# threshold sizes

n = 1024
t = 1023

print("Threshold (n, t)=", "(", n, ",", t,")")

# ring dimension, module dimension, sec dimension, tail bound

ringdim = 256

moddim = 0
secdim = 0

if balanced==1: 
	if secpar == 128: 
		moddim = 8
		secdim = 8
	if secpar == 192: 
		moddim = 11
		secdim = 11
	if secpar == 256: 
		moddim = 14
		secdim = 14

# Setting cut-off point for truncated gaussian sampling

rej_rate = 0.01

tail_bound = 0.01
tail_bound_2 = 1.01

while (2 * exp((-tail_bound**2) / 2)) > rej_rate:
    tail_bound = tail_bound + 0.01

while ((tail_bound_2**ringdim)*exp((ringdim/2)*(1-tail_bound_2^2))) > rej_rate:
    tail_bound_2 = tail_bound_2 + 0.01


# Total number of signatures

Q = 2**60

# Setting initial gaussian width of w. Set approx 38 bits. 

if balanced==1: 
	if secpar==128:
		sigma_w = 2**29
		sigma_y=2**10
	if secpar==192:
		sigma_w = 2**30
		sigma_y=2**10
	if secpar==256:
		sigma_w = 2**30
		sigma_y=2**10 


# Defining challenge space 1-norm

nu = 0
while ((2**nu) * binomial(ringdim,nu) <2**secpar):
    nu = nu + 1

# Bit-dropping parameters 

if balanced == 1:
	num = floor(log(nu,2).n()) # log2 of challspace 1-norm
	kappa_y = floor(log(sigma_w * sqrt(t * (1+secdim/moddim))/2,2)) - num
	print("kappa_y = ", kappa_y)
	kappa_w = kappa_y + num - 1
	print("kappa_w = ", kappa_w)
	print("first term bits = ",log(nu*2**(kappa_y),2).n())
	print("second term bits = ",log(2**(kappa_w+1),2).n())

# Gaussian widths for s and r, set according to Raccoon analysis



print("----------------------------------------------------------------------------------------------")

print("Assert (Rounding Correctness): nu*2^(kappa_y)+2^(kappa_w+1) <= sigma_w * sqrt[t * (1+secdim/moddim)]")
print("Assert LHS =", log(nu*2**kappa_y+2**(kappa_w+1),2).n())
print("Assert RHS =", log(sigma_w * sqrt(t * (1+secdim/moddim)),2).n())
print("Implies sigma_w can be min =", log((nu*2**kappa_y+2**(kappa_w+1))/(sqrt(t * (1+secdim/moddim))),2).n(), " bits")
sigma_w = 2**(log((nu*2**kappa_y+2**(kappa_w+1))/(sqrt(t * (1+secdim/moddim))),2)+0.0000000000001)
assert(nu*2**kappa_y+2**(kappa_w+1) <= sigma_w * sqrt(t * (1+secdim/moddim)))
print("Optimised sigma_w = ", log(sigma_w,2).n(), " bits")


print("----------------------------------------------------------------------------------------------")
print("Assert: sigma_y <= sigma_w/nu")
print("Assert LHS =", log(sigma_y,2).n())
print("Assert RHY =", log(sigma_w/nu,2).n())
assert(sigma_y <= sigma_w/nu)

# Hint-MLWE bound

B_HMLWE = Q * nu * (1+ringdim*(secpar+1+2 * log(ringdim,2))/sqrt(Q))

# Gaussian width for MLWE security

sigma = sqrt(1/(2/(n*sigma_y**2) + 2*B_HMLWE/(t*sigma_w**2)))

# 2-norm bound of signature z

B = e**(1/4) * (n * nu * sigma_y  + t * sigma_w)*sqrt(ringdim * (moddim + secdim)) + (nu * (2**kappa_y) + 2**(kappa_w+1)) * sqrt(ringdim * moddim)

#print("B_2z =", log(B,2).n())

# bound for ST-MSIS security

B_STMSIS = B + sqrt(nu) + (nu * (2**(kappa_y)) + 2**(kappa_w + 1)) * sqrt(ringdim * moddim)

#print("B_STMSIS = ", log(B_STMSIS,2).n(), "bits")

B_MSIS = B_STMSIS - nu #Following Raccoon parameter selection recommendation

#print("B_MSIS bits =", log(B_MSIS,2).n())

B_inf =  sqrt(ringdim) * tail_bound * sigma_y * sqrt(n) + tail_bound * sigma_w * sqrt(t)

# Setting signature Modulus

print("----------------------------------------------------------------------------------------------")
print("Assert (ensuring valid SIS instance): B_MSIS < (q-1)/2")
print("LHS =", log(B_MSIS,2).n())
print("Implies q can be:", log(2*B_MSIS+1,2).n(), " bits")
logq=log(2*B_MSIS+1,2)
q = round(2^logq)
while  (q % (2*ringdim) != 1) or (gcd(q-1, 2*ringdim) == 1):   #setting for NTT friendliness
	q = next_prime(q)
assert(B_MSIS < (q-1)/2)
print("q set at ", log(q,2).n(), " bits")
print("----------------------------------------------------------------------------------------------")


print("Signature parameters: ", "( N =", ringdim, ", q bits=", round(log(q,2)), ", #sigs = 2^",log(Q,2), ", chalspace 1-norm =", nu, ", kappa_y =", kappa_y, ", kappa_w =", kappa_w, ", mod dim =", moddim, ", sec dim =", secdim, ", sigma_y =", log(sigma_y,2).n(), " bits",  ", sigma_w =", log(sigma_w,2).n(), " bits", ")")


print("----------")
print("Start LWE estimation")
print("----------")

#LWE for partial signing key
y_lwe_param = LWE.Parameters(n=ringdim*secdim, q=q, Xs=ND.DiscreteGaussian(sigma_y), Xe=ND.DiscreteGaussian(sigma_y), m=ringdim*moddim, tag='Partial Signature Key LWE')

#LWE for partial commitment
w_lwe_param = LWE.Parameters(n=ringdim*secdim, q=q, Xs=ND.DiscreteGaussian(sigma_w), Xe=ND.DiscreteGaussian(sigma_w), m=ringdim*moddim, tag='Partial Commitment Key LWE')

#LWE for verification key
lwe_param = LWE.Parameters(n=ringdim*secdim, q=q, Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma), m=ringdim*moddim, tag='Combined Signature LWE')
print()

print("LWE Sec. of partial signing key")
y_lwe = LWE.estimate.rough(y_lwe_param)


print()
print("LWE Sec. of partial commitment")
w_lwe = LWE.estimate.rough(w_lwe_param)
print(w_lwe)

print()

print("LWE Sec. verification key")
r_lwe = LWE.estimate.rough(lwe_param)
print()


print("----------")
print("Start SIS estimation")
print("----------")
print()

#SIS for security
sis_param = SIS.Parameters(n=(moddim)*ringdim, q=q, length_bound=B_MSIS, m=(moddim+secdim+1)*ringdim, norm=2, tag='Signature SIS')

r_sis = SIS.estimate.rough(sis_param)

print("SIS for SelfTargetSIS")
print(r_sis)
print()


# Set parameters for encryption scheme
print("----------")
print("Start HE Parameters")
print("----------")
print()
p = q

# params for granular adjustment of security level
expansionfactor = 4
adjust_factor = 0

if secpar == 128:
	adjust_factor = 11
if secpar == 192:
	adjust_factor = 15
if secpar == 256:
	adjust_factor = 19

#secdimhat = secdim * expansionfactor
#moddimhat = moddim * expansionfactor

m = round(secdim * expansionfactor - adjust_factor)

secdimhat = m
moddimhat = m

sigma_ctx = sqrt(2/3)
sigma_B = sqrt(2/3)

B_2 = nu * n * moddimhat * tail_bound_2 * sigma_B * ringdim * tail_bound_2 * sigma_ctx + n * nu * p * tail_bound_2 * sigma_ctx*sqrt(ringdim) + p * nu * n * secdimhat * tail_bound_2 *sigma_B * ringdim * tail_bound_2 * sigma_ctx + t * tail_bound_2 * sigma_B * ringdim * tail_bound_2 * sigma_ctx * moddimhat + p * t * tail_bound_2 * sigma_ctx * sqrt(ringdim) + t * p * secdimhat * tail_bound_2 * sigma_B * ringdim *tail_bound_2 * sigma_ctx

print("B_2 bits = ", log(B_2,2).n())
print("p bits =", log(p,2).n())

sigma_tdec = sqrt(Q * B_2**2 / p**2)
print("sigma_tdec = ", sigma_tdec.n(100))

B_dec = nu * n * tail_bound_2 * sigma_y * sqrt(ringdim * secdim) + t * tail_bound_2 * sigma_w * sqrt(ringdim * secdim) + B_2 + p * t * tail_bound_2 * sigma_tdec 
q_enc = ((B_dec) * 2).n()
while  (q_enc % (2*ringdim) != 1) or (gcd(q_enc-1, 2*ringdim) == 1):
	q_enc = next_prime(q_enc)
q_enc_actual = q_enc
log_q_enc = log(q_enc_actual,2)
#print("q_enc:",N(q_enc_actual))
#print("q_enc:",str(q_enc_actual))
#print("log q_enc:", N(log_q_enc))

print("HE parameters: ", "( N =", ringdim, ", qhat bits=", N(log_q_enc,4), ", secdimhat =", secdimhat, ", moddimhat = ", moddimhat, ", m = ", m, ", secdim = ", secdim, ", moddim = ", moddim, ")")


v_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*m), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_ctx), Xe=ND.DiscreteGaussian(sigma_ctx),m=2**60, tag='HE LWE v')
u_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*secdimhat), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_ctx), Xe=ND.DiscreteGaussian(sigma_ctx),m=2**60, tag='HE LWE u')
pk_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*moddimhat), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_B), Xe=ND.DiscreteGaussian(sigma_B), tag='HE LWE Public Key')

print()
print("----------")
print("Start HE LWE estimation")
print("----------") 
print()

print("LWE Sec. of cipher text element 'v'")
v_lwe = LWE.estimate.rough(v_enc_lwe_param)
print()

print("LWE Sec. of cipher text element 'u'")
u_lwe = LWE.estimate.rough(u_enc_lwe_param)
print()

print("LWE Sec. encryption pk")
pk_lwe = LWE.estimate.rough(pk_enc_lwe_param)
print()


print("----------")
print("BDLOP Commitment Params")
print("----------")
print()

print("Computing hiding hardness:")
sig_com_rand = 0.667
A1_height = 1
A2_height = m
A1_A2_width = A1_height+A2_height+1
com_hiding_param = LWE.Parameters(n=ringdim*(A1_A2_width), q=q_enc_actual, Xs=ND.DiscreteGaussian(sig_com_rand), Xe=ND.DiscreteGaussian(sig_com_rand), m = A2_height+A1_height, tag='com hiding')
com_hiding = LWE.estimate.rough(com_hiding_param)
print(com_hiding)


print()
print("Computing binding hardness:")
B_com_rand = 16*sig_com_rand*sqrt(nu*ringdim)
print("B_com_rand =", B_com_rand.n())
sis_param = SIS.Parameters(n=(A1_height)*ringdim, q=q_enc_actual, length_bound=B_com_rand, m=ringdim*A1_A2_width, norm=2, tag='Signature SIS')
com_sis = SIS.estimate.rough(sis_param)
print(com_sis)


print()
print("----------")
print("Computing sizes")
print("----------")
print()

# Compute sizes in bits

challenge = ringdim * 1
response = ringdim * (secdim) * ceil(log(B_inf,2))
hint = ringdim * moddim * ceil(log(B/(2**kappa_w * sqrt(ringdim * (moddim))),2))

signature = challenge + response + hint
pk = ringdim * moddim * (logq - kappa_y) + secpar

print("SIG size in KB:", N(signature/8000,digits=4))
print("PK size in KB:", N(pk/8000,digits=4))
print()

##### DKG_E Passive sizing #######

h_i = 256
B_i = moddimhat * m * ringdim * log(q,2)
B_ij = moddimhat * m * ringdim * log(q,2)
S_ij = secdimhat * m * ringdim * log(q,2)
E_ij = moddimhat * m * ringdim * log(q,2)

DKG_E_passive_comcostbits = h_i + B_i + (n-1)*B_ij + (n-1)*S_ij + (n-1)*E_ij 

DKG_E_passive_comcost = (DKG_E_passive_comcostbits/8000000).n()

print("DKG_E passive com cost=", DKG_E_passive_comcost, "MB")

##### DKG_E Active sizing #######

LNP_DKGE_S=0
if secpar == 128: 
	LNP_DKGE_S= 10771152
if secpar == 192: 
	LNP_DKGE_S= 17302699
if secpar == 256: 
	LNP_DKGE_S= 24970974

pi_2_LIN_bits = m * (ringdim + (secdimhat + 2)*log(n*(2*0.675*(secdimhat+2)*ringdim*nu),2))

pi_3_LIN_bits = m * (ringdim + (secdimhat + 2)*log(t*(2*0.675*(secdimhat+2)*ringdim*nu),2))

pi_KeyGenE_bits = LNP_DKGE_S + pi_2_LIN_bits + pi_3_LIN_bits

DKG_E_allproofs_bits = pi_KeyGenE_bits
print("DKG_E proofs cost = ", (DKG_E_allproofs_bits/8000).n(), " KB")

print("DKG_E active com cost = ", ((DKG_E_allproofs_bits+DKG_E_passive_comcostbits)/8000000).n(), " MB")

##### DKG_S Passive sizing #######

h_yi = 256
yi= ringdim * secdim * log(q,2)
ctx_si = ringdim*secdimhat*log(q_enc_actual,2) + ringdim*m*log(q_enc_actual)

DKG_S_passive_comcostbits = DKG_E_passive_comcostbits + h_yi + yi + ctx_si
DKG_S_passive_comcost = (DKG_S_passive_comcostbits/8000000).n()
print("DKG_S passive com cost=", DKG_S_passive_comcost, "MB")

##### DKG_S Active sizing #######

LNP_skg_si=0
if secpar == 128: 
	LNP_skg_si= 1098331
if secpar == 192: 
	LNP_skg_si= 1238581 
if secpar == 256: 
	LNP_skg_si= 1495230

pi_LIN_skg_si_bits = 4*ringdim + log(5*(2*0.675*(secdimhat+2)*ringdim*nu),2)

print("DKG_S proofs excluding DKG_E proofs =", ((LNP_skg_si + pi_LIN_skg_si_bits)/8000).n(), " KB")

DKG_S_allproofs_bits = pi_KeyGenE_bits + LNP_skg_si + pi_LIN_skg_si_bits
print("DKG_S total proofs cost = ", (DKG_S_allproofs_bits/8000).n(), " KB")
DKG_S_active_comcostbits = DKG_S_passive_comcostbits + DKG_S_allproofs_bits
print("DKG_S total active cost = ", (DKG_S_active_comcostbits/8000000).n(), " MB")

print()

##### Passive signing sizing sizes #######

 
hash_w_size = ringdim
w_size = moddim*ringdim*logq
ctx_r_size = (secdimhat+m)*ringdim*log_q_enc
decshare_size = m*ringdim*log_q_enc

decshare_size_bits=(decshare_size/8000).n()

print("DecShare_Cost =" f"{decshare_size_bits:.2f}", "KB")

ComCost_Passive = ((round(hash_w_size + w_size + ctx_r_size + decshare_size))/8000).n()
print()
print("Passive signing communication cost = ", f"{ComCost_Passive:.2f}", "KB")
print()


###### pi_dsi proof size ######
sig_pi_dsi = 2659
com_E = N((m+1)*ringdim*log(q_enc,2)/8000,digits=4)
#print("Com(E) in KB:", com_E)
z_E = (m+3)*ringdim*log(tail_bound*sig_pi_dsi,2)
#print("z_E cost in KB: ", N(z_E/8000,digits=4))
z_S = (m+3)* m *ringdim * log(tail_bound*sig_pi_dsi,2)
#print("z_S cost in KB: ", N(z_S/8000,digits=4))
pi_LIN_dsi = N((challenge + z_E + z_S)/8000,digits=4)
#print("pi_LIN_dsi cost in KB: ", pi_LIN_dsi)


# LNP proofs calculated externally using LaZer library


LNP_ri=0
if secpar == 128: 
	LNP_ri= 123.18
if secpar == 192: 
	LNP_ri= 160.67 
if secpar == 256: 
	LNP_ri= 194.02 


LNP_dsi=0
if secpar == 128: 
	LNP_dsi=185.68
if secpar == 192: 
	LNP_dsi=238.64
if secpar == 256: 
	LNP_dsi=288.04 

pi_dsi = com_E + pi_LIN_dsi + LNP_dsi
print("Final round proof cost =: ", pi_dsi, " KB")

print("Active signing cost =", f"{ComCost_Passive + pi_dsi + LNP_ri:.2f}", "KB") 
print("Online signing cost =", f"{decshare_size_bits + pi_dsi :.2f}", "KB") 
print("Online signing cost with LaBRADOR =", f"{decshare_size_bits + com_E + pi_LIN_dsi + 66 :.2f}", "KB") 
print("Active signing cost using LaBRADOR =", f"{ComCost_Passive + 66 + LNP_ri:.2f}", "KB") 
print("Optimistic signing cost =", f"{decshare_size_bits:.2f}", "KB")




print("Done")


import math as ma
import numpy as np
import os

# Star Pasta
# Stellar Evolution Script
# Written 2021 by Nikolai Hersfeldt

# This is an implementation of the formulae in the paper "Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity", Hurley et al. 2000
# http://adsabs.harvard.edu/pdf/2000MNRAS.315..543H
# With ZAMS L and R calculations from Tout et al. 1996, as in the above paper
# https://www.semanticscholar.org/paper/Zero-age-main-sequence-radii-and-luminosities-as-of-Tout-Pols/1a2ace027e19fbba890df95360cebb5b672551e8?p2df
# And with modifications from the following papers:
# Improved NS/BH masses and EC Sne from "COMPACT REMNANT MASS FUNCTION: DEPENDENCE ON THE EXPLOSION MECHANISM AND METALLICITY", Fryer et al. 2012 (F2012)
# https://iopscience.iop.org/article/10.1088/0004-637X/749/1/91
# Additional Mass Loss mechanisms from "ON THE MAXIMUM MASS OF STELLAR BLACK HOLES", Belczynski et al. 2010 (B2010)
# https://iopscience.iop.org/article/10.1088/0004-637X/714/2/1217
# Habitable Zone Boundaries from "HABITABLE ZONES AROUND MAIN-SEQUENCE STARS: DEPENDENCE ON PLANETARY MASS", Kopparapu et al. 2014 (K2014)
# https://iopscience.iop.org/article/10.1088/2041-8205/787/2/L29/meta
# Pair-instability supernovae from "The effect of pair-instability mass loss on black-hole mergers", Belczynski et al. 2016 (B2016)


ML_on = True    #Toggles mass loss and small envelope adjustment
verbose = False #Makes script report each timestep in the command window
stop_LM = True  #Stops stars <0.7 continuing past MS
fast_SN = True  #Use remnant mass formulae from F2012 rather than old behavior; should be true per Belczynski et al. 2012, "MISSING BLACK HOLES UNVEIL THE SUPERNOVA EXPLOSION MECHANISM", https://iopscience.iop.org/article/10.1088/0004-637X/757/1/91

# PART 1: Collect user input:

print('''
Star Pasta Stellar Evolution Script
Written 2021 by Nikolai Hersfeldt
as an implementation of
"Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity", Hurley et al. 2000
 (with updates from other papers, cited in the readme)
''')

M = float(input('''Star Mass
Initial mass of the star at the start of the main sequence
High-mass stars will lose much of this over their lives
Input as solar masses (multiples of the sun's mass)
Valid between 0.08 and 300
Best accuracy between 0.8 and 150
Input: '''))

M0 = M

if M < 0.8:
    LM_prompt = input(''' Include post-main sequence in simulation?
 Accuracy not guaranteed for mass < 0.8
 and probably totally wrong for mass < 0.25
 y to include, n to stop after main sequence: ''')
    if 'y' in LM_prompt or '1' in LM_prompt:
        stop_LM = False

Z = float(input('''Star Metallicity
Portion of the star that is something other than H and He
Around 0.02 for our sun
Valid between 0.0001 and 0.03
Input: '''))


path = os.path.join(os.path.dirname(__file__), '')

##############################################################################################################

# PART 2: Determine coefficients:

def poly(z, a, b=0, c=0, d=0, e=0):
    con = a + b*z + c*z**2 + d*z**3 + e*z**4
    return con

print('Calculating coefficients...')

zeta = ma.log10(Z/0.02)
sigma = ma.log10(Z)
rho = zeta + 1

a1  = poly(zeta,  1.593890 *1000,  2.053038 *1000,  1.231226 *1000,  2.327785 *100)
a2  = poly(zeta,  2.706708 *1000,  1.483131 *1000,  5.772723 *100,   7.411230 *10)
a3  = poly(zeta,  1.466143 *100,  -1.048442 *100,  -6.795374 *10,   -1.391127 *10)
a4  = poly(zeta,  4.141960 /100,   4.564888 /100,   2.958542 /100,   5.571483 /1000)
a5  = poly(zeta,  3.426349 /10)
a6  = poly(zeta,  1.949814 *10,    1.758178,       -6.008212,       -4.470533)
a7  = poly(zeta,  4.903830)
a8  = poly(zeta,  5.212154 /100,   3.166411 /100,  -2.750074 /1000, -2.271549 /1000)
a9  = poly(zeta,  1.312179,       -3.294936 /10,    9.231860 /100,   2.610989 /100)
a10 = poly(zeta,  8.073972 /10)
a11 = poly(zeta,  1.031538,       -2.434480 /10,    7.732821,        6.460705,        1.374484)
a12 = poly(zeta,  1.043715,       -1.577474,       -5.168234,       -5.596506,       -1.299394)
a13 = poly(zeta,  7.859573 *100,  -8.542048,       -2.642511 *10,   -9.585707)
a14 = poly(zeta,  3.858911 *1000,  2.459681 *1000, -7.630093 *10,   -3.486057 *100,  -4.861703 *10)
a15 = poly(zeta,  2.888720 *100,   2.952979 *100,   1.850341 *100,   3.797254 *10)
a16 = poly(zeta,  7.196580,        5.613746 /10,    3.805871 /10,    8.398728 /100)

a11 = a11 * a14
a12 = a12 * a14

a18 = poly(zeta,  2.187715 /10,   -2.154437,       -3.768678,       -1.975518,       -3.021475 /10)
a19 = poly(zeta,  1.466440,        1.839725,        6.442199,        4.023635,        6.957529 /10)
a20 = poly(zeta,  2.652091 *10,    8.178458 *10,    1.156058 *100,   7.633811 *10,    1.950698 *10)
a21 = poly(zeta,  1.472103,       -2.947609,       -3.312828,       -9.945065 /10)
a22 = poly(zeta,  3.071048,       -5.679941,       -9.745523,       -3.594543)
a23 = poly(zeta,  2.617890,        1.019135,       -3.292551 /100,  -7.445123 /100)
a24 = poly(zeta,  1.075567 /100,   1.773287 /100,   9.610479 /1000,  1.732469 /1000)
a25 = poly(zeta,  1.476246,        1.899331,        1.195010,        3.035051 /10)
a26 = poly(zeta,  5.502535,       -6.601663 /100,   9.968707 /100,   3.599801 /100)

a17 = 10 ** max(0.097 - 0.1072 * (sigma +3), max(0.097, min(0.1461, 0.1462 + 0.1237 * (sigma + 2))))
a18 = a18 * a20
a19 = a19 * a20

a27 = poly(zeta,  9.511033 *10,    6.819618 *10,   -1.045625 *10,   -1.474939 *10)
a28 = poly(zeta,  3.113458 *10,    1.012033 *10,   -4.650511,       -2.463185)
a29 = poly(zeta,  1.413057,        4.578814 /10,   -6.850581 /100,  -5.588658 /100)
a30 = poly(zeta,  3.910862 *10,    5.196646 *10,    2.264970 *10,    2.873680)
a31 = poly(zeta,  4.597479,       -2.855179 /10,    2.709724 /10)
a32 = poly(zeta,  6.682518,        2.827718 /10,   -7.294429 /100)

a29 = a29 ** a32

a34 = poly(zeta,  1.910302 /10,    1.158624 /10,    3.348990 /100,   2.599706 /1000)
a35 = poly(zeta,  3.931056 /10,    7.277637 /100,  -1.366593 /10,   -4.508946 /100)
a36 = poly(zeta,  3.267776 /10,    1.204424 /10,    9.988332 /100,   2.455361 /100)
a37 = poly(zeta,  5.990212 /10,    5.570264 /100,   6.207626 /100,   1.777283 /100)

a33 = min(1.4, 1.5135 + 0.3769*zeta)
a33 = max(0.6355 - 0.4192*zeta, max(1.25, a33))

a38 = poly(zeta,  7.330122 /10,    5.192827 /10,    2.316416 /10,    8.346941 /1000)
a39 = poly(zeta,  1.172768,       -1.209262 /10,   -1.193023 /10,   -2.859837 /100)
a40 = poly(zeta,  3.982622 /10,   -2.296279 /10,   -2.262539 /10,   -5.219837 /100)
a41 = poly(zeta,  3.571038,       -2.223625 /100,  -2.611794 /100,  -6.359648 /1000)
a42 = poly(zeta,  1.9848,          1.1386,          3.5640   /10)
a43 = poly(zeta,  6.300    /100,   4.810    /100,   9.840    /1000)
a44 = poly(zeta,  1.200,           2.450)

a42 = min(1.25, max(1.1, a42))
a44 = min(1.3, max(0.45, a44))

a45 = poly(zeta,  2.321400 /10,    1.828075 /1000, -2.232007 /100,  -3.378734 /1000)
a46 = poly(zeta,  1.163659 /100,   3.427682 /1000,  1.421393 /1000, -3.710666 /1000)
a47 = poly(zeta,  1.048020 /100,  -1.231921 /100,  -1.686860 /100,  -4.234254 /1000)
a48 = poly(zeta,  1.555590,       -3.223927 /10,   -5.197429 /10,   -1.066441 /10)
a49 = poly(zeta,  9.7700   /100,  -2.3100   /10,   -7.5300   /100)
a50 = poly(zeta,  2.4000   /10,    1.8000   /10,    5.9500   /10)
a51 = poly(zeta,  3.3000   /10,    1.3200   /10,    2.1800   /10)
a52 = poly(zeta,  1.1064,          4.1500   /10,    1.8000   /10)
a53 = poly(zeta,  1.1900,          3.7700   /10,    1.7600   /10)

a49 = max(a49, 0.145)
a50 = min(a50, 0.306 + 0.053*zeta)
a51 = min(a51, 0.3625 + 0.062*zeta)
a52 = max(a52, 0.9)
a53 = max(a53, 1.0)
if Z > 0.01:
    a52 = min(a52, 1.0)
    a53 = min(a53, 1.1)

a54 = poly(zeta,  3.855707 /10,   -6.104166 /10,    5.676742,        1.060894 *10,    5.284014)
a55 = poly(zeta,  3.579064 /10,   -6.442936 /10,    5.494644,        1.054952 *10,    5.280991)
a56 = poly(zeta,  9.587587 /10,    8.777464 /10,    2.017321 /10)

a57 = min(1.4, 1.5135 + 0.3769*zeta)
a57 = max(0.6355 - 0.4192*zeta, max(1.25, a57))

a58 = poly(zeta,  4.907546 /10,   -1.683928 /10,   -3.108742 /10,   -7.202918 /100)
a59 = poly(zeta,  4.537070,       -4.465455,       -1.612690,       -1.623246)
a60 = poly(zeta,  1.796220,        2.814020 /10,    1.423325,        3.421036 /10)
a61 = poly(zeta,  2.256216,        3.773400 /10,    1.537867,        4.396373 /10)
a62 = poly(zeta,  8.4300   /100,  -4.7500   /100,  -3.5200   /100)
a63 = poly(zeta,  7.3600   /100,   7.4900   /100,   4.4260   /100)
a64 = poly(zeta,  1.3600   /10,    3.5200   /100)
a65 = poly(zeta,  1.564231 /1000,  1.653042 /1000, -4.439786 /1000, -4.951011 /1000, -1.216530 /1000)
a66 = poly(zeta,  1.4770,          2.9600   /10)
a67 = poly(zeta,  5.210157,       -4.143695,       -2.120870)
a68 = poly(zeta,  1.1160,          1.6600   /10)

a62 = max(0.065, a62)
if Z < 0.004:
    a63 = min(0.055, a63)
a64 = max(0.091, min(0.121, a64))
a66 = max(a66, min(1.6, -0.308 - 1.046*zeta))
a66 = max(0.8, min(0.8 - 2.0*zeta, a66))
a68 = max(0.9, min(a68, 1.0))
a68 = min(a68, a66)

if a68 > a66:
    a64 = (a58*a66**a60) / (a59 + a66**a61)

a69 = poly(zeta,  1.071489,       -1.164852 /10,   -8.623831 /100,  -1.582349 /100)
a70 = poly(zeta,  7.108492 /10,    7.935927 /10,    3.926983 /10,    3.622146 /100)
a71 = poly(zeta,  3.478514,       -2.585474 /100,  -1.512955 /100,  -2.833691 /1000)
a72 = poly(zeta,  9.132108 /10,   -1.653695 /10)
a73 = poly(zeta,  3.969331 /1000,  4.539076 /1000,  1.720906 /1000,  1.897857 /10000)
a74 = poly(zeta,  1.600,           7.640    /10,    3.322    /10)

if Z > 0.01:
    a72 = max(a72, 0.95)
a74 = max(1.4, min(a74, 1.6))

a75 = poly(zeta,  8.109    /10,   -6.282    /10)
a76 = poly(zeta,  1.192334 /100,   1.083057 /100,   1.230969,        1.551656)
a77 = poly(zeta, -1.668868 /10,    5.818123 /10,   -1.105027 *10,   -1.668070 *10)
a78 = poly(zeta,  7.615495 /10,    1.068243 /10,   -2.011333 /10,   -9.371415 /100)
a79 = poly(zeta,  9.409838,        1.522928)
a80 = poly(zeta, -2.7110   /10,   -5.7560   /10,   -8.3800   /100)
a81 = poly(zeta,  2.4930,          1.1475)

a75 = max(1.0, min(a75, 1.27))
a75 = max(a75, 0.6355 - 0.4192*zeta)
a76 = max(a76, -0.1015564 - 0.2161264*zeta - 0.05182516*zeta**2)
a77 = max(-0.3868776 - 0.5457078*zeta - 0.1463472*zeta**2, min(0.0, a77))
a78 = max(0.0, min(a78, 7.454 + 9.046*zeta))
a79 = min(a79, max(2.0, -13.3 - 18.6*zeta))
a80 = max(0.0585542, a80)
a81 = min(1.5, max(0.4, a81))

b1  = poly(zeta,  3.9700   /10,    2.8826   /10,    5.2930   /10)
b4  = poly(zeta,  9.960283 /10,    8.164393 /10,    2.383830,        2.223436,        8.638115 /10)
b5  = poly(zeta,  2.561062 /10,    7.072646 /100,  -5.444596 /100,  -5.798167 /100,  -1.349129 /100)
b6  = poly(zeta,  1.157338,        1.467883,        4.299661,        3.130500,        6.992080 /10)
b7  = poly(zeta,  4.022765 /10,    3.050010 /10,    9.962137 /10,    7.914079 /10,    1.728090 /10)

b1 = min(0.54, b1)
b2 = 10**(-4.6739 - 0.9394*sigma)
b2 = min(max(b2, -0.04167 + 55.67*Z), 0.4771 - 9329.21 * Z**2.94)
b3 = max(-0.1451, -2.2794 - 1.5175*sigma - 0.254*sigma**2)
b3 = 10**b3
if Z > 0.004:
    b3 = max(b3, 0.7307 + 14265.1 * Z**3.395)
b4 = b4 + 0.1231572*zeta**5
b6 = b6 + 0.01640687*zeta**5

b9  = poly(zeta,  2.751631 *1000,  3.557098 *100)
b10 = poly(zeta, -3.820831 /100,   5.872664 /100)
b11 = poly(zeta,  1.071738 *100,  -8.970339 *10,   -3.939739 *10)
b12 = poly(zeta,  7.348793 *100,  -1.531020 *100,  -3.793700 *10)
b13 = poly(zeta,  9.219293,       -2.005865,       -5.561309 /10)

b11 = b11**2
b13 = b13**2

b14 = poly(zeta,  2.917412,        1.575290,        5.751814 /10)
b15 = poly(zeta,  3.629118,       -9.112722 /10,    1.042291)
b16 = poly(zeta,  4.916389,        2.862149,        7.844850 /10)

b14 = b14 ** b15
b16 = b16 ** b15
b17 = 1.0
if zeta > -1.0:
    b17 = 1.0 - 0.3880523 * (zeta + 1.0) ** 2.862149

b18 = poly(zeta,  5.496045 *10,   -1.289968 *10,    6.385758)
b19 = poly(zeta,  1.832694,       -5.766608 /100,   5.696128 /100)
b20 = poly(zeta,  1.211104 *100)
b21 = poly(zeta,  2.214088 *100,   2.187113 *100,   1.170177 *10,   -2.635340 *10)
b22 = poly(zeta,  2.063983,        7.363827 /10,    2.654323 /10,   -6.140719 /100)
b23 = poly(zeta,  2.003160,        9.388871 /10,    9.656450 /10,    2.362266 /10)
b24 = poly(zeta,  1.609901 *10,    7.391573,        2.277010 *10,    8.334227)
b25 = poly(zeta,  1.747500 /10,    6.271202 /100,  -2.324229 /100,  -1.844559 /100)
b27 = poly(zeta,  2.752869,        2.729201 /100,   4.996927 /10,    2.496551 /10)
b28 = poly(zeta,  3.518506,        1.112440,       -4.556216 /10,   -2.179426 /10)

b24 = b24 ** b28
b26 = 5.0 - 0.09138012 * Z ** -0.3671407
b27 = b27 ** (2*b28)

b29 = poly(zeta,  1.626062 *100,  -1.168838 *10,   -5.498343)
b30 = poly(zeta,  3.336833 /10,   -1.458043 /10,   -2.011751 /100)
b31 = poly(zeta,  7.425137 *10,    1.790236 *10,    3.033910 *10,    1.018259 *10)
b32 = poly(zeta,  9.268325 *100,  -9.739859 *10,   -7.702152 *10,   -3.158268 *10)
b33 = poly(zeta,  2.474401,        3.892972 /10)
b34 = poly(zeta,  1.127018 *10,    1.622158,       -1.442664,       -9.474699 /10)

b31 = b31 ** b33
b34 = b34 ** b33

b36 = poly(zeta,  1.445216 /10,   -6.180219 /100,   3.093878 /100,   1.567090 /100)
b37 = poly(zeta,  1.304129,        1.395919 /10,    4.142455 /1000, -9.732503 /1000)
b38 = poly(zeta,  5.114149 /10,   -1.160850 /100)

b36 = b36 ** 4
b37 = 4.0 * b37
b38 = b38 ** 4

b39 = poly(zeta,  1.314955 *100,   2.009258 *10,   -5.143082 /10,   -1.379140)
b40 = poly(zeta,  1.823973 *10,   -3.074559,       -4.307878)
b41 = poly(zeta,  2.327037,        2.403445,        1.208407,        2.087263 /10)
b42 = poly(zeta,  1.997378,       -8.136205 /10)
b43 = poly(zeta,  1.079113 /10,    1.762409 /100,   1.096601 /100,   3.058818 /1000)
b44 = poly(zeta,  2.327409,        6.901582 /10,   -2.158431 /10,   -1.084117 /10)

b40 = max(b40, 1.0)
b41 = b41 ** b42
b44 = b44 ** 5

b46 = poly(zeta,  2.214315,       -1.975747)
b48 = poly(zeta,  5.072525,        1.146189 *10,    6.961724,        1.316965)
b49 = poly(zeta,  5.139740)

b45 = 1.0 - (2.47162*rho - 5.401682*rho**2 + 3.247361*rho**3)
if rho <= 0:
    b45 = 1.0
b47 = 1.127733*rho + 0.2344416*rho**2 - 0.3793726*rho**3

b51 = poly(zeta,  1.125124,        1.306486,        3.622359,        2.601976,        3.031270 /10)
b52 = poly(zeta,  3.349289 /10,    4.531269 /1000,  1.131793 /10,    2.300156 /10,    7.632745 /100)
b53 = poly(zeta,  1.467794,        2.798142,        9.455580,        8.963904,        3.339719)
b54 = poly(zeta,  4.658512 /10,    2.597451 /10,    9.048179 /10,    7.394505 /10,    1.607092 /10)
b55 = poly(zeta,  1.0422,          1.3156   /10,    4.5000   /100)
b56 = poly(zeta,  1.110866,        9.623856 /10,    2.735487,        2.445602,        8.826352 /10)
b57 = poly(zeta, -1.584333 /10,   -1.728865 /10,   -4.461431 /10,   -3.925259 /10,   -1.276203 /10)

b51 = b51 - 0.1343798*zeta**5
b53 = b53 + 0.4426929*zeta**5
b55 = min(0.99164 - 743.123 * Z**2.83, b55)
b56 = b56 + 0.1140142*zeta**5
b57 = b57 - 0.01308728*zeta**5

#coefficients for ZAMS luminosity and radius

lz1 = poly(zeta,   0.39704170,  -0.32913574,   0.34776688,   0.37470851,   0.09011915)
lz2 = poly(zeta,   8.52762600, -24.41225973,  56.43597107,  37.06152575,   5.45624060)
lz3 = poly(zeta,   0.00025546,  -0.00123461,  -0.00023246,   0.00045519,   0.00016176)
lz4 = poly(zeta,   5.43288900,  -8.62157806,  13.44202049,  14.51584135,   3.39793084)
lz5 = poly(zeta,   5.56357900, -10.32345224,  19.44322980,  18.97361347,   4.16903097)
lz6 = poly(zeta,   0.78866060,  -2.90870942,   6.54713531,   4.05606657,   0.53287322)
lz7 = poly(zeta,   0.00586685,  -0.01704237,   0.03872348,   0.02570041,   0.00383376)

rz1 = poly(zeta,   1.71535900,   0.62246212,  -0.92557761,  -1.16996966,  -0.30631491)
rz2 = poly(zeta,   6.59778800,  -0.42450044, -12.13339427, -10.73509484,  -2.51487077)
rz3 = poly(zeta,  10.08855000,  -7.11727086, -31.67119479, -24.24848322,  -5.33608972)
rz4 = poly(zeta,   1.01249500,   0.32699690,  -0.00923418,  -0.03876858,  -0.00412750)
rz5 = poly(zeta,   0.07490166,   0.02410413,   0.07233664,   0.03040467,   0.00197741)
rz6 = poly(zeta,   0.01077422)
rz7 = poly(zeta,   3.08223400,   0.94472050,  -2.15200882,  -2.49219496,  -0.63848738)
rz8 = poly(zeta,  17.84778000,  -7.45245690, -48.96066856, -40.05386135,  -9.09331816)
rz9 = poly(zeta,   0.00022582,  -0.00186899,   0.00388783,   0.00142402,  -0.00007671)

##############################################################################################################

# PART 3: Evolution Functions

#interpolation function; we'll use it a lot
def interp(lo_in,hi_in,lo_dep,hi_dep,var):
    res = lo_dep + (hi_dep - lo_dep) * (var - lo_in) / (hi_in - lo_in) 
    return res


# Initial mass boundaries

def f_Mhook(ze):
    Mhook = 1.0185 + 0.16015*ze + 0.0892*ze**2     #eq 1
    return Mhook

def f_MHeF(ze):
    MHeF = 1.995 + 0.25*ze + 0.087*ze**2    #eq 2
    return MHeF

def f_MFGB(z):
    MFGB = (13.048 * (z/0.02)**0.06) / (1 + 0.0012 * (0.02/z)**1.27)    #eq 3
    return MFGB

Mhook = f_Mhook(zeta)
MHeF = f_MHeF(zeta)
MFGB = f_MFGB(Z)

b46 = -1.0 * b46 * ma.log10(MHeF / MFGB)

##############################################################################################################

# Main sequence

def f_tBGB(m):
    tBGB = (a1 + a2*m**4 + a3*m**5.5 + m**7) / (a4*m**2 + a5*m**7)  #eq 4
    return tBGB

def f_tMS(m, tBGB=0, thook=0):
    if tBGB == 0:
        tBGB = f_tBGB(m)
    if thook ==0:
        thook = f_thook(m, tBGB)
    x = max(0.95, min(0.95 - 0.03 * (zeta + 0.30103), 0.99))  #eq 6
    tMS = max(thook, x*tBGB)    #eq 5
    return tMS

def f_thook(m, tBGB=0):
    if tBGB == 0:
        tBGB = f_tBGB(m)
    mu = max(0.5, 1.0 - 0.01 * max(a6/m**a7, a8 + a9/m**a10))   #eq 7
    thook = mu * tBGB
    return thook

#LZAMS and RZAMS from Tout et al. 1996

def f_LZAMS(m):
    LZAMS = (lz1*m**5.5 + lz2*m**11) / (lz3 + m**3 + lz4*m**5 + lz5*m**7 + lz6*m**8 + lz7*m**9.5)
    return LZAMS

def f_RZAMS(m):
    RZAMS = (rz1*m**2.5 + rz2*m**6.5 + rz3*m**11 + rz4*m**19 + rz5*m**19.5) / (rz6 + rz7*m**2 + rz8*m**8.5 + m**18.5 + rz9*m**19.5)
    return RZAMS

#Back to Hurley et al(2000)

def f_LTMS(m):
    LTMS = (a11*m**3 + a12*m**4 + a13*m**(a16+1.8)) / (a14 + a15*m**5 + m**a16)     #eq 8
    return LTMS

def f_RTMS_1(m):
    RTMS = (a18 + a19*m**a21) / (a20 + m**a22)      #eq 9a
    return RTMS

def f_RTMS_2(m):
    c1 = -8.672073 /100
    RTMS = (c1*m**3 + a23*m**a26 + a24*m**(a26+1.5)) / (a25 + m**5)     #eq 9b
    return RTMS

def f_RTMS(m, RZAMS=0):
    if m <= a17:
        RTMS = f_RTMS_1(m)
        if m < 0.5:
            if RZAMS == 0:
                RZAMS = f_RZAMS(m)
            RTMS = max(RTMS, 1.5 * RZAMS)
    elif m >= a17 + 0.1:
        RTMS = f_RTMS_2(m)
    else:
        lo = f_RTMS_1(a17)
        hi = f_RTMS_2(a17 + 0.1)
        RTMS = interp(a17, a17+0.1, lo, hi, m)
    return RTMS

def f_LBGB(m):
    c2 = 9.301992
    c3 = 4.637345
    LBGB = (a27*m**a31 + a28*m**c2) / (a29 + a30*m**c3 + m**a32)    #eq 10
    return LBGB


def main_seq(m, t):
    tBGB = f_tBGB(m)
    thook = f_thook(m, tBGB)
    tMS = f_tMS(m, tBGB, thook)
    tau = t / tMS   #eq 11
    
    LZAMS = f_LZAMS(m)
    LTMS = f_LTMS(m)
    RZAMS = f_RZAMS(m)
    RTMS = f_RTMS(m)
    
    tau1 = min(1.0, t/thook)    #eq 14
    tau2 = max(0.0, min(1.0, (t - (1.0 - 0.01) * thook) / (0.01 * thook)))      #eq 15
    
    if m <= Mhook:      #eq 16
        dL = 0.0
    elif m < a33:
        B = min(a34 / a33**a35, a36 / a33**a37)
        dL = B * ((m - Mhook) / (a33 - Mhook))**0.4
    else:
        dL = min(a34 / m**a35, a36 / m**a37)
    
    if m <= Mhook:      #eq 17
        dR = 0.0
    elif m <= a42:
        dR = a43 * ((m - Mhook) / (a42 - Mhook))**0.5
    elif m < 2.0:
        B = (a38 + a39*2**3.5) / (a40*2**3 + 2**a41) - 1.0
        dR = a43 + (B - a43) * ((m - a42) / (2.0 - a42))**a44
    else:
        dR = (a38 + a39*m**3.5) / (a40*m**3 + m**a41) - 1.0
    
    if Z > 0.0009 or m <= 1.0:      #eq 18
        eta = 10
    else:
        if m >= 1.1:
            eta = 20
        else:
            eta = interp(1.0, 1.1, 10, 20, m)
    
    if m < 0.5:     #eq 19b
        aL = a49
    elif m < 0.7:
        aL = a49 + 5.0 * (0.3 - a49) * (m - 0.5)
    elif m < a52:
        aL = 0.3 + (a50 - 0.3) * (m - 0.7) / (a52 - 0.7)
    elif m < a53:
        aL = a50 + (a51 - a50) * (m - a52) / (a53 - a52)
    elif m < 2.0:
        B = (a45 + a46*2**a48) / (2**0.4 + a27*2**1.9)
        aL = a51 + (B - a51) * (m - a53) / (2.0 - a53)
    else:
        aL = (a45 + a46*m**a48) / (m**0.4 + a27*m**1.9)     #eq 19a
    
    BL = max(0.0, a54 - a55*m**a56)     #eq 20
    if m > 0.0 and BL > 0.0:
        B = a54 - a55*a57**a56
        BL = max(0.0, B - 10.0 * (m - a57) * B)

    if m < 0.5:     #eq 21b
        aR = a62
    elif m < 0.65:
        aR = a62 + (a63 - a62) * (m - 0.5) / 0.15
    elif m < a68:
        aR = a63 + (a64 - a63) * (m - 0.65) / (a68 - 0.65)
    elif m < a66:
        B = (a58*a66**a60) / (a59 + a66**a61)
        aR = a64 + (B - a64) * (m - a68) / (a66 - a68)
    elif m <= a67:
        aR = (a58*m**a60) / (a59 + m**a61)      #eq 21a
    else:
        C = (a58*a67**a60) / (a59 + a67**a61)
        aR = C + a65 * (m - a67)
    
    if m <= 1.0:    #eq 22b
        BR = 1.06
    elif m < a74:
        BR = 1.06 + (a72 - 1.06) * (m - 1.0) / (a74 - 1.06)
    elif m < 2.0:
        B = (a69*2**3.5) / (a70 + 2**a71)
        BR = a72 + (B - a72) * (m - a74) / (2.0 - a74)
    elif m <= 16:
        BR = (a69*m**3.5) / (a70 + m**a71)      #eq 22a
    else:
        C = (a69*16**3.5) / (a70 + 16**a71)
        BR = C + a73 * (m - 16.0)
    BR = BR - 1
    
    if m <= 1.0:    #eq 23
        y = a76 + a77 * (m - a78)**a79
    elif m <= a75:
        B = a76 + a77 * (1 - a78)**a79
        y = B + (a80 - B) * ((m - 1.0) / (a75 - 1.0))**a81
    elif m < a75 + 0.1:
        if a75 == 1.0:
            C = a76 + a77 * (1 - a78)**a79
        else:
            C = a80
        y = C - 10.0 * (m - a75) * C
    else:
        y = 0.0
    
    LMS_1 = aL*tau + BL*tau**eta + (ma.log10(LTMS / LZAMS) - aL - BL) * tau**2 - dL*(tau1**2 - tau2**2)     #eq 12
    LMS = LZAMS * 10**LMS_1
    
    RMS_1 = aR*tau + BR*tau**10 + y*tau**40 + (ma.log10(RTMS / RZAMS) - aR - BR - y) * tau**3 - dR*(tau1**3 - tau2**2)      #eq 13
    RMS = RZAMS * 10**RMS_1
    if m < 0.1:
        X = 0.76 - 3.0*Z
        RMS = max(RMS, 0.0258 * (1.0 + X)**(5/3) * m**(-1/3))   #eq 24
    
    return LMS, RMS

##############################################################################################################

#Hertzsprung Gap

def f_LEHG(m):
    if m < MFGB:
        LEHG = f_LBGB(m)
    else:
        LEHG = f_LHeI(m)
    return LEHG

def f_REHG(m):
    if m < MFGB:
        LBGB = f_LBGB(m)
        REHG = f_RGB(m, LBGB)
    else:
        REHG = f_RHeI(m)
    return REHG

def f_McEHG(m):
    if m < MHeF:    #eq 28
        LBGB = f_LBGB(m)
        McEHG = f_McGB(m, LBGB)
    elif m < MFGB:
        McEHG = f_McBGB_IM(m)
    else:
        McEHG = f_McHeI(m)
    return McEHG

def f_McHG(m, tau=0, t=0):
    if tau==0:
        tBGB = f_tBGB(m)
        tMS = f_tMS(m, tBGB)
        tau = (t - tMS) / (tBGB - tMS)      #eq 25
    McEHG = f_McEHG(m)
    rho_1 = (1.586 + m**5.25) / (2.434 + 1.02 * m**5.25)    #eq 29
    McHG = ((1-tau)*rho_1 + tau) * McEHG    #eq 30
    return McHG

def hertz_gap(m, mt, t):
    tBGB = f_tBGB(m)
    tMS = f_tMS(m, tBGB)
    tau = (t - tMS) / (tBGB - tMS)      #eq 25
    
    LEHG = f_LEHG(m)
    LTMS = f_LTMS(m)
    LHG = LTMS * (LEHG / LTMS)**tau     #eq 26
    
    REHG = f_REHG(mt)
    RTMS = f_RTMS(mt)
    RHG = RTMS * (REHG / RTMS)**tau     #eq 27

    McHG = f_McHG(m, tau, t)
    
    return LHG, RHG, McHG

##############################################################################################################

# First Giant Branch

def f_LGB(m, Mc, p=0, q=0, B=0, D=0):
    if p == 0:
        p = f_p(m)
    if q == 0:
        q = f_q(m)
    if B == 0:
        B = f_B(m)
    if D == 0:
        D = f_D(m)
    L = min(B * Mc**q, D * Mc**p)   #eq 37
    return L

def f_McGB(m, L, p=0, q=0, B=0, D=0):
    if p == 0:
        p = f_p(m)
    if q == 0:
        q = f_q(m)
    if B == 0:
        B = f_B(m)
    if D == 0:
        D = f_D(m)
    Mx = f_Mx(m, p, q, B, D)
    Lx = f_LGB(m, Mx, p, q, B, D)
    if L < Lx:
        Mc = (L / B)**(1/q)     #eq 37 reversed
    else:
        Mc = (L / D)**(1/p)
    return Mc

def f_Mx(m, p=0, q=0, B=0, D=0):
    if p == 0:
        p = f_p(m)
    if q == 0:
        q = f_q(m)
    if B == 0:
        B = f_B(m)
    if D == 0:
        D = f_D(m)
    Mx = (B / D)**(1 / (p - q))     #eq 38
    return Mx

def f_p(m):
    if m <= MHeF:
        p = 6
    elif m >= 2.5:
        p = 5
    else:
        p = interp(MHeF, 2.5, 6, 5, m)
    return p

def f_q(m):
    if m <= MHeF:
        q = 3
    elif m >= 2.5:
        q = 2
    else:
        q = interp(MHeF, 2.5, 3, 2, m)
    return q

def f_B(m):
    B = max(30000.0,500.0 + 17500.0*m**0.6)
    return B

def f_D(m):
    D0 = 5.37 + 0.135*zeta
    if m <= MHeF:
        D_1 = D0
    elif m >= 2.5:
        D_1 = max(-1.0, 0.975*D0 - 0.18*m, 0.5*D0 - 0.06*m)
    else:
        hi = max(-1.0, 0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)
        D_1 = interp(MHeF, 2.5, D0, hi, m)
    D = 10**D_1
    return D

def f_McGB_1(t, tinf1, tx, tinf2, p, q, B, D, AH):
    if t <= tx:
        McGB_1 = ((p - 1) * AH * D * (tinf1 - t))**(1/(1-p))    #eq 39
    else:
        McGB_1 = ((q - 1) * AH * B * (tinf2 - t))**(1/(1-q))
    return McGB_1

def f_tinf1(m, tBGB=0, LBGB=0, p=0, D=0, AH=0):
    if tBGB == 0:
        tBGB = f_tBGB(m)
    if LBGB == 0:
        LBGB = f_LBGB(m)
    if p == 0:
        p = f_p(m)
    if D == 0:
        D = f_D(m)
    if AH == 0:
        AH = f_AH(m)
    tinf1 = tBGB + 1/((p-1)*AH*D) * (D / LBGB)**((p-1)/p)   #eq 40
    return tinf1

def f_tx(m, tinf1=0, Lx=0, tBGB=0, LBGB=0, p=0):
    if tBGB == 0:
        tBGB = f_tBGB(m)
    if LBGB == 0:
        LBGB = f_LBGB(m)
    if p == 0:
        p = f_p(m)
    if tinf1 == 0:
        tinf1 = f_tinf1(m, tBGB, LBGB, p)
    if Lx == 0:
        Mx = f_Mx(m, p)
        Lx = f_LGB(m, Mx, p)
    tx = tinf1 - (tinf1 - tBGB) * (LBGB / Lx)**((p-1)/p)    #eq 41
    return tx

def f_tinf2(m, tx=0, Lx=0, q=0, B=0, AH=0):
    if q == 0:
        q = f_q(m)
    if B == 0:
        B = f_B(m)
    if AH == 0:
        AH = f_AH(m)
    if tx == 0:
        tx = f_tx(m)
    if Lx == 0:
        p = f_p(m)
        Mx = f_Mx(m, p)
        Lx = f_LGB(m, Mx, p)
    tinf2 = tx + 1 / ((q - 1) * AH * B) * (B / Lx)**((q-1)/q)   #eq 42
    return tinf2

def f_tHeI(m, tinf1=0, tinf2=0, p=0, q=0, B=0, D=0, AH=0, LHeI=0, Lx=0):
    if LHeI == 0:
        LHeI = f_LHeI(m)
    if p == 0:
        p = f_p(m)
    if q == 0:
        q = f_q(m)
    if B == 0:
        B = f_B(m)
    if D == 0:
        D = f_D(m)
    if AH == 0:
        AH = f_AH(m)
    if Lx == 0:
        Mx = f_Mx(m, p, q, B, D)
        Lx = f_LGB(m, Mx, p, q, B, D)
    if LHeI <= Lx:
        if tinf1 == 0:
            tinf1 = f_tinf1(m)
            tHeI = tinf1 - 1/((p-1)*AH*D) * (D / LHeI)**((p-1)/p)   #eq 43
        if tinf2 == 0:
            tinf2 = f_tinf2(m)
            tHeI = tinf1 - 1/((q-1)*AH*B) * (B / LHeI)**((q-1)/q)
    return tHeI

def f_AH(m):
    AH_1 = max(-4.8, min(-5.7 + 0.8*m, -4.1 + 0.14*m))
    AH = 10**AH_1
    return AH

def f_McBGB_IM(m):
    c1 = 9.20925 / 100000
    c2 = 5.402216
    LBGB_MHeF = f_LBGB(MHeF)
    Mc = f_McGB(m, LBGB_MHeF)
    C = Mc**4 - c1*MHeF**c2
    McBAGB = f_McBAGB(m)
    McBGB = min(0.95 * McBAGB, (C + c1*m**c2)**0.25)    #eq 44
    return McBGB

def f_McHeI(m, LHeI=0):
    if m < MHeF:
        if LHeI == 0:
            LHeI = f_LHeI(m)
        McHeI = f_McGB(m, LHeI)
    else:
        c1 = 9.20925 / 100000
        c2 = 5.402216
        LHeI_MHeF = f_LHeI(MHeF)
        Mc = f_McGB(m, LHeI_MHeF)
        C = Mc**4 - c1*MHeF**c2
        McBAGB = f_McBAGB(m)
        McHeI = min(0.95 * McBAGB, (C + c1*m**c2)**0.25)
    return McHeI

def f_RGB(m, L):
    A = min(b4*m**-b5, b6*m**-b7)
    RGB = A * (L**b1 + b2*L**b3)    #eq 46
    return RGB

def giant_branch(m, mt, t):
    tBGB = f_tBGB(m)
    
    tHeI = f_tHeI(m)
    
    tau = (t - tBGB) / (tHeI - tBGB)
    
    p = f_p(m)
    q = f_q(m)
    B = f_B(m)
    D = f_D(m)
    
    AH = f_AH(m)
    
    LBGB = f_LBGB(m)
    
    Mx = f_Mx(m, p, q, B, D)
    Lx = f_LGB(m, Mx, p, q, B, D)
    
    tinf1 = f_tinf1(m, tBGB, LBGB, p, D, AH)
    tx = f_tx(m, tinf1, Lx, tBGB, LBGB, p)
    tinf2 = f_tinf2(m, tx, Lx, q, B, AH)
    
    McGB_1 = f_McGB_1(t, tinf1, tx, tinf2, p, q, B, D, AH)
    
    LGB = f_LGB(m, McGB_1, p, q, B, D)
    
    if m < MHeF:
        McGB = McGB_1
    else:
        McBGB = f_McBGB_IM(m)
        McHeI = f_McHeI(m)
        McGB = McBGB + (McHeI - McBGB)*tau      #eq 45
    
    RGB = f_RGB(mt, LGB)
    
    return LGB, RGB, McGB

##############################################################################################################

# Core Helium Burning

def f_LHeI(m):
    if m < MHeF:
        LHeI_MHeF = (b11 + b12*MHeF**3.8) / (b13 + MHeF**2)
        alpha = (b9*MHeF**b10 - LHeI_MHeF) / LHeI_MHeF
        LHeI = b9*m**b10 / (1 + alpha * ma.exp(15 * (m - MHeF)))    #eq 49
    else:
        LHeI = (b11 + b12*m**3.8) / (b13 + m**2)
    return LHeI

def f_RHeI(m, LHeI=0, RmHe=0):
    if m <= MFGB:
        if LHeI == 0:
            LHeI = f_LHeI(m)
        RHeI = f_RGB(m, LHeI)
    elif m >= max(MFGB, 12.0):
        if RmHe == 0:
            RmHe = f_RmHe(m)
        RHeI = RmHe
    else:
        if LHeI == 0:
            LHeI = f_LHeI(m)
        if RmHe == 0:
            RmHe = f_RmHe(m)
        mu_1 = ma.log10(m / 12.0) / ma.log10(MFGB / 12.0)
        RHeI = RmHe * (f_RGB(m, LHeI) / RmHe)**mu_1     #eq 50
    return RHeI

def f_LminHe(m, LHeI=0):
    if LHeI == 0:
        LHeI = f_LHeI(m)
    c = b17 / MFGB**0.1 + (b16*b17 - b14) / MFGB**(b15+0.1)
    LminHe = LHeI * (b14 + c*m**(b15+0.1)) / (b16 + m**b15)     #eq 51
    return LminHe

def f_mu(m, Mc):
    mu = (m - Mc) / (MHeF - Mc)     #eq 52
    return mu

def f_LZAHB(m, Mc, mu=0):
    if mu == 0:
        mu = f_mu(m, Mc)
    LZHe_Mc = f_LZHe(Mc)
    LminHe_MHeF = f_LminHe(MHeF)
    alpha2 = (b18 + LZHe_Mc - LminHe_MHeF) / (LminHe_MHeF - LZHe_Mc)
    LZAHB = LZHe_Mc + (1 + b20) / (1 + b20*mu**1.6479) * b18*mu**b19 / (1 + alpha2 * ma.exp(15*(m - MHeF)))    #eq 53
    return LZAHB

def f_RZAHB(m, Mc, mu=0, LZAHB=0):
    if mu == 0:
        mu = f_mu(m, Mc)
    if LZAHB == 0:
        LZAHB = f_LZAHB(m, Mc, mu)
    f = ((1.0 + b21) * mu**b22) / (1.0 + b21*mu**b23)
    RZAHB = (1-f) * f_RZHe(Mc) + f * f_RGB(m, LZAHB)    #eq 54
    return RZAHB

def f_RmHe(m, LZAHB=0):
    if m >= MHeF:
        RmHe = (b24*m + (b25*m)**b26 * m**b28) / (b27 + m**b28)     #eq 55
    else:
        assert LZAHB != 0, 'No LZAHB passed to RmHe calculation'
        RmHe_MHeF = (b24*MHeF + (b25*MHeF)**b26 * MHeF**b28) / (b27 + MHeF**b28)
        LZAHB_MHeF = f_LminHe(MHeF)     #LZAHB and LminHe should be continuous at MHeF
        RmHe = f_RGB(m, LZAHB) * (RmHe_MHeF / (f_RGB(MHeF, LZAHB_MHeF)))**(m/MHeF)
    return RmHe

def f_LBAGB(m):
    if m < MHeF:
        LBAGB_MHeF = (b31 + b32*MHeF**(b33+1.8)) / (b34 + MHeF**b33)
        alpha3 = (b29*MHeF**b30 - LBAGB_MHeF) / LBAGB_MHeF
        LBAGB = b29*m**b30 / (1 + alpha3 * ma.exp(15 * (m - MHeF)))    #eq 56
    else:
        LBAGB = (b31 + b32*m**(b33+1.8)) / (b34 + m**b33)
    return LBAGB

def f_tHe(m, tBGB=0, Mc=0, mu=0):
    if tBGB == 0:
        tBGB = f_tBGB(m)
    if m < MHeF:
        assert Mc != 0, 'No Mc passed to tHe calculation'
        if mu == 0:
            mu = f_mu(m, Mc)
        tHeMS_Mc = f_tHeMS(Mc)
        alpha4 = (tBGB * (b41*MHeF**b42 + b43*MHeF**5) / (b44 + MHeF**5) - b39) / b39
        tHe = (b39 + (tHeMS_Mc - b39) * (1 - mu)**b40) * (1 + alpha4 * ma.exp(15 * (m - MHeF)))     #eq 57
    else:
        tHe = tBGB * (b41*m**b42 + b43*m**5) / (b44 + m**5)
    return tHe

def f_taubl(m):
    if m < MHeF:
        taubl = 1
    elif m <= MFGB:
        alphabl = (1 - b45 * (MHeF / MFGB)**0.414) * (ma.log10(MHeF / MFGB))**-b46
        taubl = b45 * (m / MFGB)**0.414 + alphabl * (ma.log10(m / MFGB))**b46       #eq 58
        taubl = taubl.real      #It doesn't quite seem right that there are complex numbers in here, but I can't see how else you'd interpret the formulae, and the imaginary component does end up being negligible
    else:
        taubl = (1 - b47) * f_fbl(m) / f_fbl(MFGB)
    if m > MFGB and f_RmHe(m) >= f_RAGB(m, f_LHeI(m)):    #Not explicitly stated in the paper but implied as necessary
        taubl = 0
    taubl = max(0, min(1, taubl))
    return taubl

def f_fbl(m, RmHe=0, LHeI=0):
    if RmHe == 0:
        RmHe = f_RmHe(m)
    if LHeI == 0:
        LHeI = f_LHeI(m)
    fbl = m**b48 * (1 - RmHe / f_RAGB(m, LHeI))**b49
    return fbl

def f_McBAGB(m):
    McBAGB = (b36*m**b37 + b38)**0.25   #eq 66
    return McBAGB

def core_he_burn(m, mt, t, Mc):
    tHeI = f_tHeI(m)
    tBGB = f_tBGB(m)
    
    LHeI = f_LHeI(m)
    
    McHeI = f_McHeI(m, LHeI)
    
    McBAGB = f_McBAGB(m)
    
    loops = 0
    McCHeB = 0
    while abs(McCHeB - Mc) / Mc > 0.999:    #McCHeB and tau are interdependent, so we iterate calculating them until there's <0.1% variance in the McCHeB calculation between iterations
        McCHeB = Mc
        
        mu = f_mu(m, McCHeB)
        tHe = f_tHe(m, tBGB, McCHeB, mu)
        tau = (t - tHeI) / tHe
        
        Mc = (1 - tau) * McHeI + tau * McBAGB   #eq 67
        loops += 1
        assert loops < 1000, 'Too many loops in McCHeB calculation'     #prevents infinite loops
    McCHeB = Mc
    taubl = f_taubl(m)
    if m < MHeF or m > MFGB:
        taux = 0
    else:
        taux = 1 - taubl
        
    if m < MFGB:    #I think this check should prevent some errors, but I'm not sure honestly
        LZAHB = f_LZAHB(m, McCHeB, mu)
    
    LminHe = f_LminHe(m, LHeI)
    
    if m < MHeF:
        Lx = LZAHB      #eq 59
    elif m < MFGB:
        Lx = LminHe
    else:
        Lx = LHeI
    
    if m < MFGB:
        RmHe = f_RmHe(mt, LZAHB)
    else:
        if mt < MHeF:   #this feels a bit hacky but it prevents errors on high-mass stars
            LZAHB = f_LZAHB(mt, McCHeB, mu)
            RmHe = f_RmHe(mt, LZAHB)
        else:
            RmHe = f_RmHe(mt)
    RHeI = f_RHeI(mt, LHeI, RmHe)
    
    if mt < MHeF:
        RZAHB = f_RZAHB(mt, McCHeB, mu, LZAHB)
        Rx = RZAHB      #eq 60
    elif mt < MFGB:
        Rx = f_RGB(mt, LminHe)
    else:
        Rx = RHeI
    
    xi = min(2.5, max(0.4, RmHe / Rx))
    lamb = ((tau - taux) / (1 - taux))**xi  #eq 62
    if taux > tau:
        lamb_1 = ((taux - tau) / taux)**3      #eq 63
    
    LBAGB = f_LBAGB(m)
    
    if taux <= tau:
        LCHeB = Lx * (LBAGB / Lx)**lamb     #eq 61
    else:
        LCHeB = Lx * (LHeI / Lx)**lamb_1
    
    Rmin = min(RmHe, Rx)
    
    if m > MFGB:
        tauy = taubl
    else:
        tauy = 1
    
    if m <= MFGB:
        Ly = LBAGB
    else:
        if taux <= tau:
            Ly = Lx * (LBAGB / Lx)**lamb
        else:
            Ly = Lx * (LHeI / Lx)**lamb_1
    
    Ry = f_RAGB(mt, Ly)
    
    if tau < taux:
        RCHeB = f_RGB(mt, LCHeB)    #eq 64
    elif tau > tauy:
        RCHeB = f_RAGB(mt, LCHeB)
    elif tau == 0.0 and tauy == 0.0:     #prevents div0 errors for stars too large to have a blue loop
        RCHeB = RHeI
    else:
        rho_1 = (ma.log(Ry / Rmin))**(1/3) * ((tau - taux) / (tauy - taux)) - (ma.log(Rx / Rmin))**(1/3) * ((tauy - tau) / (tauy - taux))   #eq 65
        RCHeB = Rmin * ma.exp(abs(rho_1)**3)
    
    return LCHeB, RCHeB, McCHeB

##############################################################################################################

# Asymptotic Giant Branch

def f_McDU(McBAGB):
    if McBAGB <= 0.8:
        McDU = McBAGB
    else:
        McDU = 0.44 * McBAGB + 0.448
    return McDU

def f_RAGB(m, L):
    if m >= MHeF:
        b50 = b55 * b3
        A = min(b51*m**-b52, b53*m**-b54)
    elif m <= MHeF - 0.2:
        b50 = b3
        A = b56 + b57*m
    else:
        b50_lo = b3
        b50_hi = b55 * b3
        b50 = interp(MHeF-0.2, MHeF, b50_lo, b50_hi, m)
        A_lo = b56 + b57 * (MHeF - 0.2)
        A_hi = min(b51*MHeF**-b52, b53*MHeF**-b54)
        A = interp(MHeF-0.2, MHeF, A_lo, A_hi, m)
    RAGB = A * (L**b1 + b2*L**b50)      #eq 74
    return RAGB

def Asymptotic(m, mt, t):
    p = f_p(m)
    q = f_q(m)
    B = f_B(m)
    D = f_D(m)
    AH = f_AH(m)
    AHe = 7.66 / 100000     #eq 68
    
    Mx = f_Mx(m, p, q, B, D)
    Lx = f_LGB(m, Mx, p, q, B, D)
    
    McBAGB = f_McBAGB(m)
    
    tBGB = f_tBGB(m)
    tHeI = f_tHeI(m)
    tHe = f_tHe(m, tBGB, McBAGB)
    
    tBAGB = tHeI + tHe
    
    LBAGB = f_LBAGB(m)
    
    McDU = f_McDU(McBAGB)
    LDU = f_LGB(m, McDU, p, q, B, D)
    
    tinf1 = f_tinf1(m, tBAGB, LBAGB, p, D, AHe)
    tx = f_tx(m, tinf1, Lx, tBAGB, LBAGB, p)
    tinf2 = f_tinf2(m, tx, Lx, q, B, AHe)
    
    if LDU <= Lx:
        tDU = tinf1 - 1 / ((p - 1) * AHe * D) * (D / LDU)**((p-1)/p)    #eq 70
    else:
        tDU = tinf2 - 1 / ((q - 1) * AHe * B) * (B / LDU)**((q-1)/q)
    
    if t <= tDU or McBAGB > 2.25:
        late = False
        McHe = McBAGB
        McCO = f_McGB_1(t, tinf1, tx, tinf2, p, q, B, D, AHe)
        LAGB = f_LGB(m, McCO, p, q, B, D)
    else:
        late = True     #signals that the TPAGB has begun
        AHHe = (AH * AHe) / (AH + AHe)      #eq 71
        if LDU <= Lx:
            tinf1 = f_tinf1(m, tDU, LDU, p, D, AHHe)
            tx = f_tx(m, tinf1, Lx, tDU, LDU, p)
            tinf2 = f_tinf2(m, tx, Lx, q, B, AHHe)
            
            McCO_1 = f_McGB_1(t, tinf1, tx, tinf2, p, q, B, D, AHHe)
        else:
            tinf2_1 = tDU + 1 / ((q - 1) * AHHe * B) * (B / LDU)**((q-1)/q)     #eq 72
            McCO_1 = ((q - 1) * AHHe * B * (tinf2_1 - t))**(1/(1-q))
        LAGB = f_LGB(m, McCO_1, p, q, B, D)
        
        lamb = min(0.9, 0.3 + 0.001*m**5)   #eq 73
        
        McCO = McDU + (1 - lamb) * (McCO_1 - McDU)
        McHe = McCO
    RAGB = f_RAGB(mt, LAGB)
        
    return LAGB, RAGB, McHe, McCO, late

##############################################################################################################

#Naked Helium Stars

def f_LZHe(m):
    LZHe = 15262*m**10.25 / (m**9 + 29.54*m**7.5 + 31.18*m**6 + 0.0469)     #eq 77
    return LZHe

def f_RZHe(m):
    RZHe = 0.2391*m**4.6 / (m**4 + 0.162*m**3 + 0.0065)     #eq 78
    return RZHe

def f_tHeMS(m):
    tHeMS = (0.4129 + 18.81*m**4 + 1.853*m**6) / m**6.5     #eq 79
    return tHeMS

def f_LHeMS(m, tau, LZHe=0):
    if LZHe == 0:
        LZHe = f_LZHe(m)
    alpha = max(0, 0.85 - 0.08*m)   #eq 82
    LHeMS = LZHe * (1 + 0.45*tau + alpha*tau**2)    #eq 80
    return LHeMS

def f_RHeMS(mt, tau, RZHe=0):
    if RZHe == 0:
        RZHe = f_RZHe(mt)
    beta = max(0, 0.4 - 0.22 * ma.log10(mt))    #eq 83
    RHeMS = RZHe * (1 + beta*tau - beta*tau**6)     #eq 81
    return RHeMS

def f_Mcmax(m):
    Mcmax = min(1.45 * m - 0.31, m)
    return Mcmax

def f_LHeGB(m, McCO):
    B_1 = 4.1 * 10000
    D_1 = 5.5 * 10000 / (1 + 0.4*m**4)
    
    LHeGB = min(B_1 * McCO**3, D_1 * McCO**5)   #eq 84
    return LHeGB

def f_RHeGB(mt, LHeGB, LTHe=0, RZHe=0):
    if RZHe == 0:
        RZHe = f_RZHe(mt)
    if LTHe == 0:
        LTHe = f_LHeMS(mt, 1)
    lamb = 500 * (2 + mt**5) / mt**2.5      #eq 87
    R1 = RZHe * (LHeGB /  LTHe)**0.2 + 0.02 * (ma.exp(LHeGB / lamb) - ma.exp(LTHe / lamb))  #eq 86
    R2 = 0.08 * LHeGB**0.75     #eq 88
    
    RHeGB = min(R1, R2)     #eq 85
    
    if R1 < R2:
        late = True     #signals that the HeGB has begun
    else:
        late = False
    return RHeGB, late
    

def he_main_seq(m, mt, t):
    LZHe = f_LZHe(m)
    RZHe = f_RZHe(mt)
    tHeMS = f_tHeMS(m)
    
    tau = t / tHeMS
    
    LHeMS = f_LHeMS(m, tau, LZHe)
    
    RHeMS = f_RHeMS(mt, tau, RZHe)
    
    return LHeMS, RHeMS

def he_giant_branch(m, mt, t):
    p = f_p(m)
    q = f_q(m)
    B = f_B(m)
    D = f_D(m)
    
    AHe = 7.66 / 100000
    
    LTHe = f_LHeMS(m, 1)
    
    Mx = f_Mx(m, p, q, B, D)
    Lx = f_LGB(m, Mx, p, q, B, D)
    
    tHeMS = f_tHeMS(m)
    
    tinf1 = f_tinf1(m, tHeMS, LTHe, p, D, AHe)
    tx = f_tx(m, tinf1, Lx, tHeMS, LTHe, p)
    tinf2 = f_tinf2(m, tx, Lx, q, B, AHe)
    
    McCO = f_McGB_1(t, tinf1, tx, tinf2, p, q, B, D, AHe)
    
    LHeGB = f_LHeGB(m, McCO)
    
    RZHe = f_RZHe(mt)
    RHeGB, late = f_RHeGB(mt, LHeGB, LTHe, RZHe)
    
    return LHeGB, RHeGB, McCO

##############################################################################################################

# Remnants

def f_McSN(McBAGB):
    McSN = max(1.44, 0.773 * McBAGB - 0.35)     #eq 75
    return McSN

def f_LWD(m, t, A):
    LWD = 635 * m * Z**0.4 / (A * (t + 0.1))**1.4   #eq 90
    return LWD

def f_RWD(m):
    MCh = 1.44
    RNS = 1.4 / 100000
    
    RWD = max(RNS, 0.0115 * ma.sqrt((MCh / m)**(2/3) - (m / MCh)**(2/3)))   #eq 91
    return RWD

def white_dwarf(m, t, A):
    LWD = f_LWD(m, t, A)
    
    RWD = f_RWD(m)
    
    return LWD, RWD

def f_MNS(McSN):    #Superseded by additions from F2012
    MNS = 1.17 + 0.09 * McSN
    return MNS

def neutron(m, t):
    LNS = 0.02 * m**(2/3) / max(t,0.1)**2   #eq 93
    RNS = 1.4 / 100000
    
    return LNS, RNS

def black_hole(m,t):
    LBH = 10**-10     #eq 96
    RBH = 4.24 / 1000000 * m    #eq 94
    
    return LBH, RBH

##############################################################################################################

# PART 4: Mass Loss and Envelope Compensation

# Small-Envelope Behavior

def small_env(m, m0, Mc, McCO, L, R, stage, t):
    if ML_on == False or stage < 2 or stage == 7 or stage > 9:
        L_1 = L
        R_1 = R
        Rcr = 0
    else:
        if stage > 7:
            Mcmax = f_Mcmax(m0)
            mu = 5 * ((Mcmax - McCO) / Mcmax)     #eq 98
        else:
            L0 = 70000
            kappa = -0.5
            mu = ((m - Mc) / m) * min(5.0, max(1.2, (L / L0)**kappa))   #eq 97
        if stage < 4:
            if m >= MHeF:    # I'm pretty sure the original paper has this the wrong way round due to a typo based on the surrounding text
                Lc = f_LZHe(Mc)
                Rc = f_RZHe(Mc)
            else:
                Lc = f_LWD(Mc, 0, 4)
                Rc = f_RWD(min(Mc,1.44))    #I've added the restriction of only passing masses under the Chandrasekhar limit, to prevent errors
        elif stage == 4:
            tHeI = f_tHeI(m0)
            tBGB = f_tBGB(m0)
            tHe = f_tHe(m0, tBGB, Mc)
            tau = (t - tHeI) / tHe
            Lc = f_LHeMS(Mc, tau)
            Rc = f_RHeMS(Mc, tau)
        elif stage == 5:
            Lc = f_LHeGB(Mc, McCO)
            Rc, NA = f_RHeGB(Mc, Lc)
        else:
            Lc = f_LWD(McCO, 0, 15)
            Rc = f_RWD(min(McCO,1.44))
        Rc = min(R*0.9999, Rc)
        if stage == 4 or stage == 5:
            Rcr = Rc
        else:
            Rcr = 5 * Rc
        if mu < 1.0:
            b = 0.002 * max(1, 2.5 / m)     #eq 103
            c = 0.006 * max(1, 2.5 / m)     #eq 104
            q = ma.log(R / Rc)      #eq 105
            s = ((1 + b**3) * (mu / b)**3) / (1 + (mu / b)**3)      #eq 101
            r = ((1 + c**3) * (mu / c)**3 * mu**(0.1/q)) / (1 + (mu / c)**3)    #eq 102
            L_1 = Lc * (L / Lc)**s      #eq 99
            R_1 = Rc * (R / Rc)**r      #eq 100
        else:
            L_1 = L
            R_1 = R
            
    return L_1, R_1, Rcr

# Mass Loss

def mass_loss(m, Mc, McCO, L, R, stage):
    if ML_on == False or stage > 9 or stage == 0:
        ML = 0.0
    else:
        MR = 0.0
        MVW = 0.0
        MNJ = 0.0
        MWR = 0.0
        MOB = 0.0
        if stage > 2:
            eta = 0.5
            MR = eta * (4/10**13) * eta * L * R / m     #eq 106
        if stage == 5 or stage == 6:
            P0 = 10**(min(3.3, -2.07 - 0.9 * ma.log10(m) + 1.94 * ma.log10(R)))
            MVW = 10**(-11.4 + 0.0125 * (P0 - 100 * max(m - 2.5, 0.0)))
            MVW = min(MVW, (1.36 / 10**9) * L)
        if L > 4000:
            MNJ = (9.6/10**15) * (Z/0.02)**0.5 * R**0.81 * L**1.24 * m**0.16
        if stage > 1:
            if stage > 6:
                mu = 0.0
            else:
                L0 = 70000
                kappa = -0.5
                if stage == 6:
                    Mc = McCO
                mu = ((m - Mc) / m) * min(5.0, max(1.2, (L / L0)**kappa))   #eq 97
            if mu < 1.0:
                #MWR = 10**-13 * L**1.5 * (1.0 - mu)
                MWR = 10**-13 * L**1.5 * (Z/0.02)**0.86 * (1.0 - mu)    #eq 9 from B2010
        if stage < 7:   #Additional OB mass-loss per B2010
            Teff = 5778 * ma.sqrt(ma.sqrt(L) / R)
            if Teff >= 12500 and Teff <= 50000:
                if Teff < 25000:    #eq 6 from B2010
                    MOB = 10**(-6.688 + 2.210 * ma.log10(L/100000) - 1.339 * ma.log10(m/30) - 1.601 * ma.log10(1.3/2) + 0.85 * zeta + 1.07 * ma.log10(Teff/20000))
                else:   #eq 7 from B2010
                    MOB = 10**(-6.697 + 2.194 * ma.log10(L/100000) - 1.313 * ma.log10(m/30) - 1.226 * ma.log10(2.6/2) + 0.85 * zeta + 0.933 * ma.log10(Teff/40000) - 10.92 * ma.log10(Teff/40000)**2)
        ML = max (MR, MVW, MNJ, MWR, MOB)
        if L > 600000 and (R * L)**0.5 / 100000 > 1.0:
            if stage < 7:
                ML = 1.5/10000  #eq 8 from B2010
            else:
                MLBV = 0.1 * ((R * L)**0.5 / 100000 - 1.0)**3 * (L / 600000 - 1.0)
                ML = ML + MLBV
        if stage > 6:
            ML = max(MR, MWR)
        ML = ML * 1e6
    return ML

##############################################################################################################

# PART 5: Control Functions

#determines the timestep to use in the next simulation step
def timestep(m, ML, t, stage, Mc=0, McCO=0, mt=0):
    if stage == 1:
        tMS = f_tMS(m)
        dtk = tMS / 100
        dte = tMS - t
    elif stage == 2:
        if m > MFGB:
            tBGB = f_tHeI(m)
        else:
            tBGB = f_tBGB(m)
        tMS = f_tMS(m, tBGB)
        dtk = (tBGB - tMS) / 20
        dte = tBGB - t
    elif stage == 3:
        tBGB = f_tBGB(m)
        p = f_p(m)
        q = f_q(m)
        B = f_B(m)
        D = f_D(m)
        AH = f_AH(m)
        LBGB = f_LBGB(m)
        Mx = f_Mx(m, p, q, B, D)
        Lx = f_LGB(m, Mx, p, q, B, D)
        tinf1 = f_tinf1(m, tBGB, LBGB, p, D, AH)
        tx = f_tx(m, tinf1, Lx, tBGB, LBGB, p)
        tinf2 = f_tinf2(m, tx, Lx, q, B, AH)
        if t <= tx:
            dtk = (tinf1 - t) / 50
        else:
            tinf2 = f_tinf2(m, tx)
            dtk = (tinf2 - t) / 50
        tHeI = f_tHeI(m)
        dte = tHeI - t
    elif stage == 4:
        assert Mc != 0, 'No Mc passed to timestep function in CHeB stage'
        tBGB = f_tBGB(m)
        tHe = f_tHe(m, tBGB, Mc)
        dtk = tHe / 50
        tHeI = f_tHeI(m)
        dte = tHeI + tHe - t
    elif stage == 5 or stage == 6:
        p = f_p(m)
        q = f_q(m)
        B = f_B(m)
        D = f_D(m)
        AHe = 7.66 / 100000
        Mx = f_Mx(m, p, q, B, D)
        Lx = f_LGB(m, Mx, p, q, B, D)
        McBAGB = f_McBAGB(m)
        tBGB = f_tBGB(m)
        tHeI = f_tHeI(m)
        tHe = f_tHe(m, tBGB, McBAGB)
        tBAGB = tHeI + tHe
        LBAGB = f_LBAGB(m)
        tinf1 = f_tinf1(m, tBAGB, LBAGB, p, D, AHe)
        tx = f_tx(m, tinf1, Lx, tBAGB, LBAGB, p)
        tinf2 = f_tinf2(m, tx, Lx, q, B, AHe)
        McDU = f_McDU(McBAGB)
        LDU = f_LGB(m, McDU, p, q, B, D)
        if LDU <= Lx:
            tDU = tinf1 - 1 / ((p - 1) * AHe * D) * (D / LDU)**((p-1)/p)    #eq 70
        else:
            tDU = tinf2 - 1 / ((q - 1) * AHe * B) * (B / LDU)**((q-1)/q)
        if t < tDU:
            dte = tDU - t
        else:
            dte = 10**8 - t
            AH = f_AH(m)
            AHHe = (AH * AHe) / (AH + AHe)      #eq 71
            if LDU <= Lx:
                tinf1 = f_tinf1(m, tDU, LDU, p, D, AHe)
                tx = f_tx(m, tinf1, Lx, tDU, LDU, p)
                tinf2 = f_tinf2(m, tx, Lx, q, B, AHe)
            else:
                tinf2 = tDU + 1 / ((q - 1) * AHHe * B) * (B / LDU)**((q-1)/q)
        if t <= tx:
            dtk = (tinf1 - t) / 50
        else:
            dtk = max(1e-6,abs(tinf2 - t) / 50)   #Feels a bit wrong that t can be over tinf2--it usually doesn't cause issues but I do need to limit the minimum timestep
        if stage == 6:
            dtk = min(dtk, 5e-3)
    elif stage == 7:
        tHeMS = f_tHeMS(m)
        dtk = tHeMS / 20
        dte = tHeMS - t
    elif stage == 8 or stage == 9:
        p = f_p(m)
        q = f_q(m)
        B = f_B(m)
        D = f_D(m)
        AHe = 7.66 / 100000
        LTHe = f_LHeMS(m, 1)
        Mx = f_Mx(m, p, q, B, D)
        Lx = f_LGB(m, Mx, p, q, B, D)
        tHeMS = f_tHeMS(m)
        tinf1 = f_tinf1(m, tHeMS, LTHe, p, D, AHe)
        tx = f_tx(m, tinf1, Lx, tHeMS, LTHe, p)
        if t <= tx:
            dtk = (tinf1 - t) / 50
        else:
            tinf2 = f_tinf2(m, tx, Lx, q, B, AHe)
            dtk = max(1e-6,abs(tinf2 - t) / 50)
        dte = 10**8 - t
    else:
        dtk = max(0.1, 10 * t)
        dte = 10**8 - t
    
    if ML > 0:
        dtml = 0.01 * m/ML      #limits mass change due to mass loss to 1% per timestep
        if stage > 9:
            dtn = 10**8
        elif stage > 6:
            dtn = (mt - McCO) / ML
        else:
            dtn = (mt - Mc) / ML
    else:
        dtml = 10**8
        dtn = 10**8
    #dtk = max(1e-3,dtk)     #previous versions needed this to keep the TPAGB from getting too long, but it seems okay now
    dt = max(0.0, min(dtk,dte,dtml,dtn))
    return dt


#controls evolution of star between stages
def evolve(m, t, stin, mt=0, Mc=0, McCO=0, late=False):
    m0 = m
    t1 = t
    if stin == 8 or stin == 9:
        McSN = f_McSN(m)
        McBAGB = m
    else:
        McBAGB = f_McBAGB(m)
        McSN = f_McSN(McBAGB)
    if stin > 6:
        Mc = McCO
    if stin > 9:
        Mc = 0.0
        McCO = 0.0
        McSN = 1.0
    if McBAGB >= 1.83 and McBAGB <= 2.25 and McCO > 1.38: #Per F2012; the text is not totally clear on the implementation so this is a best guess
        print(' Electron-Capture Supernova!')
        t1 = 0.0
        m0 = 1.26078    #Solution for eq 13 in F2012 with Mrembar of 1.38
        stout = 13
        print(' Neutron Star')
    elif McCO >= McSN:    #Supernova or collapse event; includes several modifications from F2012
        if stin == 8 or stin == 9:
            McBAGB = m
        else:
            McBAGB = f_McBAGB(m)
        t1 = 0.0
        if McBAGB < 1.83:   #Altered from 1.6 per F2012
            stout = 15
            print(' Supernova!')
            print(' No Remnant')
            m0 = 0
        else:
            direct = False
            if fast_SN:
                Mproto = 1      #eq 15 from F2012
                if McSN < 2.5:  #eq 16 from F2012
                    Mfb = 0.2
                elif McSN < 6:
                    Mfb = 0.286 * McSN - 0.514
                elif McSN < 7:
                    Mfb = m - 1
                    direct = True
                elif McSN < 11:
                    ka1 = 0.25 - 1.275 / (m - 1)
                    kb1 = -11 * ka1 + 1
                    Mfb = (m - 1) * (ka1 * McSN + kb1)
                else:
                    Mfb = m - 1
                    direct = True
            else:
                if McSN < 4.82:     #eq 10 from F2012
                    Mproto = 1.50
                elif McSN < 6.31:
                    Mproto = 2.11
                elif McSN < 6.75:
                    Mproto = 0.69 * McSN - 2.26
                else:
                    Mproto = 0.37 * McSN - 0.07
                if McSN < 5:        #eq 11 from F2012
                    Mfb = 0
                elif McSN < 7.6:
                    Mfb = (m - Mproto) * (0.378 * McSN - 1.889)
                else:
                    Mfb = m - Mproto
                    direct = True
            Mrembar = Mproto + Mfb  #eq 12 from F2012
            if Mrembar > 2.96875:   #corresponding to a maximum NS mass of 2.5
                if stin > 6:
                    MHe = m
                else:
                    MHe = Mc
                if MHe > 45 and MHe < 135:  #pair-instability mechanisms per B2016
                    if MHe > 65:
                        m0 = 0
                        stout = 15
                        print(' Pair-Instability Supernova!')
                        print(' No Remnant')
                    else:
                        m0 = 40.5
                        stout = 14
                        print(' Pair-Instability Pulsation Supernova!')
                        print(' Black Hole')
                else:
                    if direct:
                        print(' Direct Collapse')
                    else:
                        print(' Supernova!')
                    stout = 14
                    m0 = 0.9 * Mrembar  #eq 14 from F2012
                    print(' Black Hole')
            else:
                print(' Supernova!')
                stout = 13
                m0 = (ma.sqrt(1 + 0.3 * Mrembar) - 1) / 0.15    #solution for eq 13 from F2012
                print(' Neutron Star')
    elif Mc >= mt:      #Transitions to White Dwarf or Naked Helium star
        m0 = Mc
        t1 = 0.0
        if stin == 6:
            McBAGB = f_McBAGB(m)
            if McBAGB < 1.83:   #Altered from 1.6 per F2012
                stout = 11
                print(' C/O White Dwarf')
            else:
                stout = 12
                print(' O/Ne White Dwarf')
        elif stin == 5:
            stout = 8
            Mcmax = f_Mcmax(Mc)
            if McCO >= Mcmax:   #Some low-mass, high-metallicity stars can reach this point with a core mass already above the expected point of He fusion cessation
                stout = 11      #The paper does not bring up the possibility of stars skipping straight from EAGB to WD but I can't see how else to interpret this
                print(' C/O White Dwarf')
            else:
                print(' Naked Helium HG')
                p = f_p(Mc)
                q = f_q(Mc)
                B = f_B(Mc)
                D = f_D(Mc)
                AHe = 7.66 / 100000
                Mx = f_Mx(m, p, q, B, D)
                LTHe = f_LHeMS(Mc, 1)
                tHeMS = f_tHeMS(Mc)
                tinf1 = f_tinf1(Mc, tHeMS, LTHe, p, D, AHe)
                if Mc > Mx:
                    Lx = f_LGB(Mc, Mx, p, q, B, D)
                    tx = f_tx(m, tinf1, Lx, tHeMS, LTHe, p)
                    tinf2 = f_tinf2(Mc, tx, Lx, q, B, AHe)
                    t1 = tinf2 - McCO**(1-q) / ((q - 1) * AHe * B)      #reverse of eq 39
                else:
                    t1 = tinf1 - McCO**(1-p) / ((p - 1) * AHe * D)
        elif stin == 4:
            Mcmax = f_Mcmax(Mc)
            if McCO >= Mcmax:   #As above
                stout = 11
                print(' C/O White Dwarf')
            else:
                stout = 7
                print(' Naked Helium MS')
                tHeI = f_tHeI(m)
                tBGB = f_tBGB(m)
                tHe = f_tHe(m, tBGB, Mc)
                tHeMS = f_tHeMS(Mc)
                t1 = ((t - tHeI) / tHe) * tHeMS     #eq 76
        else:
            if m < MHeF:
                stout = 10
                print(' He White Dwarf')
            else:
                stout = 7
                print(' Naked Helium MS')
    elif stin == 1:
        tMS = f_tMS(m)
        tMS1 = f_tMS(mt)
        t1 = t * tMS1 / tMS
        m0 = mt
        if t >= tMS:
            if stop_LM and m < 0.8:
                stout = 0   #I've opted to use "0" to designate the cessation of modelling for low-mass stars, rather than their main sequence as in Hurley et al. 2000
                print(' Main Sequence Concluded; Post-MS not simulated')
            else:
                stout = 2
                print(' Hertzsprung Gap')
        else:
            stout = 1
    elif stin == 2:
        if m0 > MFGB:
            tBGB = f_tHeI(m0)
        else:
            tBGB = f_tBGB(m0)
        Mc1 = f_McHG(mt, 0, t)  #Only reduce initial mass if it would not cause an unphysical decrease in core mass
        if Mc1 < Mc:
            tBGB1 = tBGB
        else:
            m0 = mt
            if mt > MFGB:
                tBGB1 = f_tHeI(mt)
            else:
                tBGB1 = f_tBGB(mt)
            tMS = f_tMS(m)
            tMS1 = f_tMS(mt)
            t1 = tMS1 + (t - tMS) * (tBGB1 - tMS1) / (tBGB - tMS)
        if t1 >= tBGB1:
            if m0 > MFGB:
                stout = 4
                print(' Core Helium Burning')
            else:
                stout = 3
                print(' First Giant Branch')
        else:
            stout = 2
    elif stin == 3:
        tHeI = f_tHeI(m)
        if t >= tHeI:
            stout = 4
            print(' Core Helium Burning')
            if m < MHeF:
                m0 = mt
                tHeI = f_tHeI(m0)
                t1 = tHeI
        else:
            stout = 3
    elif stin == 4:
        tHeI = f_tHeI(m)
        tBGB = f_tBGB(m)
        tHe = f_tHe(m, tBGB, Mc)
        if t >= tHeI + tHe:
            stout = 5
            print(' Early Asymptotic Giant Branch')
        else:
            stout = 4
    elif stin == 5:
        if late:
            stout = 6
            print(' Thermally Pulsating Asymptotic Giant Branch')
        else:
            stout = 5
    elif stin == 7:
        Mcmax = f_Mcmax(m)
        if McCO >= Mcmax:
            stout = 11
            t1 = 0.0
            print(' C/O White Dwarf')
        tHeMS = f_tHeMS(m)
        tHeMS1 = f_tHeMS(mt)
        t1 = t * tHeMS1 / tHeMS
        if t1 >= tHeMS1:
            stout = 8
            print(' Naked Helium HG')
        else:
            stout = 7
            m0 = mt
    elif stin == 8 or stin == 9:
        Mcmax = f_Mcmax(m)
        if McCO >= Mcmax:
            stout = 11
            t1 = 0.0
            print(' C/O White Dwarf')
        elif late:
            stout = 9
            if stin == 8:
                print(' Naked Helium GB')
        else:
            stout = 8
    else:
        stout = stin
    return stout, m0, t1


#Checks that change over the last timestep hasn't been too rapid
def retry_check(mt, m0, Mc, McCO, R, R1, stage, stagei):
    good = 1
    if abs(R1 - R) > 0.1 * R and stage == stagei:   #limits radius change to 10% of current radius in any timestep, except during a stage change
        good = 0                                    #to save a bit on calculation time, the compared radii are both pre-small envelope adjustment
    if stage < 10:
        if stage == 8 or stage == 9:
            McSN = f_McSN(m0)
        else:
            McBAGB = f_McBAGB(m0)
            McSN = f_McSN(McBAGB)
        if McCO > McSN * 1.01:     #prevents overshooting supernova
            good = 0
        if mt < Mc * 0.99:     #prevents overshooting envelope loss
            good = 0
        elif mt < McCO * 0.99:
            good = 0
    return good
    

#Computes extra data and stores it all to numpy array
def data_store(datain, step, stage, t, mt, Mc, McCO, ML, L, R, Rc):
    data = np.append(datain, [[step, stage, t, mt, Mc, McCO, ML, L, R, Rc]], 0)
    return data

#General formula for HZ calculation
def f_HZ(L, THZ, S, ka, kb, kc, kd):
    Te = THZ - 5780
    Seff = S + ka*Te + kb*Te**2 + kc*Te**3 + kd*Te**4  #eq 4 from K2014
    dist = (L / Seff)**0.5  #eq 5 from K2015
    return dist
    

#Computes extra data for output
def data_add(L, R):
    Teff = 5778 * ma.sqrt(ma.sqrt(L) / R)
    THZ = max(2600, min(7200, Teff))    #Fits from K2014 are only given for 2600-7200 K, so the Seff values at these bounds are crudely extrapolated outwards
    hzoptin = f_HZ(L, THZ, 1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15)   #Table 1 from K2014
    hzconin = f_HZ(L, THZ, 1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15)
    hzconout = f_HZ(L, THZ, 0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16)
    hzoptout = f_HZ(L, THZ, 0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16)
    return Teff, hzoptin, hzconin, hzconout, hzoptout
    
    

#Simulates one timestep
def sim_step(m0, mti, ML, Mci, McCOi, R1, t0, dt, stagei, late):
    good = 0
    while good == 0:
        t1 = t0 + dt
        mt = mti - ML * dt
        stage, m, t1 = evolve(m0, t1, stagei, mt, Mci, McCOi, late)
        late = False
        if m < mt:
            mt = m
        if stage == 1:
            L, R = main_seq(m, t1)
            Mc = 0.0
            McCO = 0.0
        elif stage == 2:
            L, R, Mc = hertz_gap(m, mt, t1)
            McCO = 0.0
        elif stage == 3:
            L, R, Mc = giant_branch(m, mt, t1)
            McCO = 0.0
        elif stage == 4:
            L, R, Mc = core_he_burn(m, mt, t1, Mci)
            McCO = 0.0
        elif stage == 5 or stage == 6:
            L, R, Mc, McCO, late = Asymptotic(m, mt, t1)
        elif stage == 7:
            L, R = he_main_seq(m, mt, t1)
            Mc = mt
            McCO = 0.0
        elif stage == 8 or stage == 9:
            L, R, McCO = he_giant_branch(m, mt, t1)
            Mc = mt
        elif stage > 9 and stage < 13:
            if stage == 10:
                A = 4
            elif stage == 11:
                A = 15
            else:
                A = 17
            L, R = white_dwarf(m, t1, A)
            Mc = mt
            if stage > 10:
                McCO = 0
            else:
                McCO = mt
        elif stage == 13:
            L, R = neutron(m, t1)
            Mc = mt
            McCO = mt
        elif stage == 14:
            L, R = black_hole(m, t1)
            Mc = mt
            McCO = mt
        else:
            L = 0
            R = 0
            Mc = 0
            McCO = 0
        good = retry_check(mti-ML*dt, m0, Mci, McCOi, R, R1, stage, stagei)
        if good == 0:
            if dt < 1e-6:
                good = 1
            else:
                dt = dt / 2
        if Mc > mt:     #prevents star magically gaining mass when core mass overshoots current mass
            Mc = mt
        if McCO > mt:
            McCO = mt
        if stage > 6 and stage < 9:
            Mcmax = f_Mcmax(m)
            if McCO > Mcmax:
                McCO = Mcmax
    return m, mt, Mc, McCO, t1, dt, L, R, stage, late
            
#Core simulation loop
def sim_run(m):
    print('Evolving Star...')
    stage = 1
    late = False
    step = 0
    t = 0.0
    t1 = t
    dt = 0.0
    m0 = M
    mt = m0
    ML = 0.0
    Mc = 0.0
    McCO = 0.0
    R = f_RZAMS(m0)
    R1 = R
    stagei = 0
    data = np.empty([0,15])
    print(' Main Sequence')
    while t < 10**8:
        stagei = stage
        m0, mt, Mc, McCO, t1, dt, L, R1, stage, late = sim_step(m0, mt, ML, Mc, McCO, R1, t1, dt, stage, late)
        t = t + dt
        if stage != stagei:
            print('  (' + str(t) + ' myr)')
        if verbose == True:
            print(t)
        L, R, Rcr = small_env(mt, m0, Mc, McCO, L, R1, stage, t1)
        ML = mass_loss(mt, Mc, McCO, L, R, stage)
        if stage == 0 or stage == 15:
            Teff, hzoptin, hzconin, hzconout, hzoptout = 0, 0, 0, 0, 0
        else:
            Teff, hzoptin, hzconin, hzconout, hzoptout = data_add(L, R)
        data = np.append(data, [[step, stage, t, mt, Mc, McCO, ML, L, R, Rcr, Teff, hzoptin, hzconin, hzconout, hzoptout]], 0)
        dt = timestep(m0, ML, t1, stage, Mc, McCO, mt)
        step += 1
        if step == 2001:
            print(' WARNING: very long output')
            print('  you may need to extend the lists in the "Organized Data" tab of the starpasta_out spreadsheet')
        elif step == 100001:
            print(' WARNING: extremely long output')
            print('  may be beyond length that the starpasta_out spreadsheet can handle')
        if stage == 0 or stage == 15:
            break
    print('Simulation Complete')
    print('Writing data...')
    savename = path + 'Z' + str(Z) + '_M' + str(m) + '.csv'
    np.savetxt(savename, data, delimiter=",")
    print('Output saved to ' + str(savename))
    return
    
    
        
        
        
        
        


##############################################################################################################

# PART 6: Main Routine

sim_run(M)
input('')




















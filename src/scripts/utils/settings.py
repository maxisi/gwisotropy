# minimum component mass allowed [MSUN]
MMIN = 3.
MMAX = 150.
ZMIN = 1.
ZMAX = 3.

# maximum false alarm threshold (1/yr)
FARMAX_YR = 1.0
# maximum false alarm threshold (1/s)
one_per_year = 1/(365.25*24.0*3600.0)
FARMAX_S = FARMAX_YR * one_per_year

# number of samples to draw per event
NSAMP = 2048

# standard deviation of hyperprior
SIGMA_RAW = 0.4

# number of isotropized validation runs
NITER_ISO = 8

NSTART_SEL = 4
NITER_SEL = 12

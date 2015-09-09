from math import *
import numpy as np

import PASTIS_v6
import PASTIS_v6.models.RV as rvmodule

FWHM1 = 10
FWHM2 = 12
fluxratio = 1e-2
deltarv = 10

nrepet = 100

# Oversampled CCF parameters
stepccf = 0.001
rvccf = np.arange(-250., 250., stepccf)

sigma1 = FWHM1/(2 * sqrt(2*log(2)))
sigma2 = FWHM2/(2 * sqrt(2*log(2)))

# Prepare array for fit results
dt = [
    ("vspan", np.float64),
    ("wspan", np.float64),
    ("fwhm", np.float64),
    ]

# Deal with non-uniform naming
fitfunc = {'vspan' : rvmodule.vspan,
           'bis' : rvmodule.fit_BIS,
           
           
           
results = np.empty(nrepet, dtype=dt)

for i in range(nrepet):
    # Introduce random shift in sampling
    x = rvccf + (np.random.rand() - 0.5)

    # Produced blended line profile
    ccf1 = 1 - 0.3 * np.exp(x**2 / (2 * sigma1**2))
    ccf2 = flxratio * (1 - 0.3 * np.exp((x - deltarv)**2 / (2 * sigma2**2)))
    ccf = (ccf1 + ccf2) / (1 + fluxratio)

    # Bin on fixed rv points

    # Fit
    
    vspan
    wspan
    BiGauss
    RV
    FWHM

from math import *
import numpy as np
import pylab as p
import time

import PASTIS_v6
import PASTIS_v6.tools
import PASTIS_v6.models.RV as rvmodule

FWHM1 = 10
FWHM2 = 12
fluxratio = 1e-2
deltarv = 10

nrepet = 100
binsize = np.logspace(0, 4, 10).astype(np.int)

# Value to use in jiggling the RV sampling [in km/s].
jiggle = 1.0

# Oversampled CCF parameters
stepccf = 0.001
rvccf = np.arange(-50., 50., stepccf)

sigma1 = FWHM1/(2 * sqrt(2*log(2)))
sigma2 = FWHM2/(2 * sqrt(2*log(2)))

# Prepare array for fit results
dt = [
    #("rv", np.float64),
    #("fwhm", np.float64),
    ("vspan", np.float64),
    ("bis", np.float64),
    ("bigauss", np.float64),
    ]

results = np.empty((len(binsize), nrepet), dtype=dt)

# Computing times initialized to zero
comptime = np.zeros(len(binsize), dtype=dt)
bintime = np.zeros(len(binsize), dtype=np.float64)

# Deal with non-uniform naming
fitfunc = {#'rv' : rvmodule.fitgauss
           'vspan' : rvmodule.vspan,
           'bis' : rvmodule.fit_BIS,
           'bigauss' : rvmodule.fit_BiGauss,
           }

# Array of shifted rvs
x = rvccf + (jiggle*np.random.rand(nrepet).reshape(nrepet, 1) - 0.5)

# Produced blended line profiles
ccf1 = 1 - 0.3 * np.exp(-x**2 / (2 * sigma1**2))
ccf2 = fluxratio * (1 - 0.3 * np.exp(-(x - deltarv)**2 /
                                     (2 * sigma2**2)))
ccf = (ccf1 + ccf2) / (1 + fluxratio)

for j, binfactor in enumerate(binsize):

    print('Binning CCF to {:.3f} km/s '
          'resolution.'.format(binfactor*stepccf))

    for i in range(len(x)):

        # Bin on fixed rv points
        ti = time.time()
        xb = PASTIS_v6.tools.rebin(x[i], binfactor)
        ccfb = PASTIS_v6.tools.rebin(ccf[i], binfactor)
        bintime[j] += time.time() - ti
        
        # Produce educated starting point for fits.
        p0 = [0.4, xb[np.argmin(ccfb)], 10.0]

        # Fit for each diagnostic
        for diag in results.dtype.names:
            ti = time.time()
            results[diag][j, i] = fitfunc[diag](xb, ccfb, p0)
            # Add time to counter
            comptime[diag][j] += time.time() - ti

xx = binsize*stepccf

f1 = p.figure()
ax = f1.add_subplot(111)

f2 = p.figure()
axr = f2.add_subplot(111)

f3 = p.figure()
axt = f3.add_subplot(111)

colordict = {'vspan' : 'r',
             'bis' : 'b',
             'bigauss' : 'k'}

for diag in results.dtype.names:
    ax.loglog(xx, results[diag].std(axis=1), color=colordict[diag], lw=2,
              label=diag)
    axr.loglog(xx, results[diag].std(axis=1)/np.median(results[diag][0]),
               color=colordict[diag], lw=2, label=diag)

    # Benchmark
    y = np.log10((comptime[diag] + bintime)/nrepet)
    par = np.polyfit(np.log10(xx)[:-2], y[:-2], 1)

    axt.loglog(xx, comptime[diag]/nrepet, color=colordict[diag], lw=1,
               ls=':')

    axt.loglog(xx, (comptime[diag] + bintime)/nrepet,
               color=colordict[diag], lw=2, ls='-',
               label='{0} (O {1:.2f})'.format(diag, par[0]))
    
axt.loglog(xx, bintime/nrepet, color='0.55', lw=2, ls='-', label='binning')
    
for aa in (ax, axr, axt):
    aa.legend(loc=0)
    aa.set_xlabel('Sampling size [km/s]', fontsize=16)

ax.set_ylabel('Absolute precision [km/s]', fontsize=16)
axr.set_ylabel('Relative precision', fontsize=16)
axt.set_ylabel('Computing time', fontsize=16)

p.show()

import numpy as np
from astropy.modeling import fitting, models

# generate data
y, x = np.indices((51, 51))
data = 3 * np.exp(-((x - 24) ** 2) / (2 * 3**2) - (y - 26) ** 2 / (2 * 4**2))

# guess parameters
y0, x0 = np.unravel_index(np.argmax(data), data.shape)
sigmax = np.mean(np.std(data, axis=0))
sigmay = np.mean(np.std(data, axis=1))
amp = np.max(data)
# initialise 2D gaussian
w = models.Gaussian2D(amp, x0, y0, sigmax, sigmay)
# run fit
fit_w = fitting.LMLSQFitter()
yi, xi = np.indices(data.shape)
g = fit_w(w, xi, yi, data)
# results
print(w)
# predictions
predictions = g(xi, yi)

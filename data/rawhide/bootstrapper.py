## bootstrapper.py ############################################################
## tools to create a bootstrap for a linear regression model, and plot this  ##
## in matplotlib ###############################################################
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def plotBootstrapsOnDataPlot(pllt, x, y, strapColor='grey', regressColor='red'):
	
	
	
	# Extend x data to contain another row vector of 1s
	x = np.asarray(x)
	y = np.asarray(y)

	X = np.vstack([x, np.ones(len(x))]).T

	plt.figure(figsize=(12,8))
	for i in range(0, 500):
		sample_index = np.random.choice(range(0, len(y)), len(y))

		X_samples = X[sample_index]
		y_samples = y[sample_index]    

		lr = LinearRegression()
		lr.fit(X_samples, y_samples)
		plt.plot(x, lr.predict(X), color=strapColor, alpha=0.1, zorder=1)

	lr = LinearRegression()
	lr.fit(X, y)
	plt.plot(x, lr.predict(X), color=regressColor, zorder=5)	

if(__name__ == "__main__"):
	## Create toy data and bootstrap it to verify everything is working
	## correctly
	
	x = np.linspace(0, 10, 20)
	y = x + (np.random.rand(len(x)) * 10)
	plotBootstrapsOnDataPlot(plt, x, y)
	plt.scatter(x, y, marker='+', color='blue', zorder=5)

	plt.savefig('boostrapDemo.png')

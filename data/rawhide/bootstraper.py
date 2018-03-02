## bootstrapper.py ############################################################
###############################################################################

import numpy as no
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression



if(__name__ == "__main__"):
	# Create toy data 
	x = np.linspace(0, 10, 20)
	y = x + (np.random.rand(len(x)) * 10)

	# Extend x data to contain another row vector of 1s
	X = np.vstack([x, np.ones(len(x))]).T

	plt.figure(figsize=(12,8))
	for i in range(0, 500):
		sample_index = np.random.choice(range(0, len(y)), len(y))

		X_samples = X[sample_index]
		y_samples = y[sample_index]    

		lr = LinearRegression()
		lr.fit(X_samples, y_samples)
		plt.plot(x, lr.predict(X), color='grey', alpha=0.2, zorder=1)

	plt.scatter(x,y, marker='o', color='orange', zorder=4)

	lr = LinearRegression()
	lr.fit(X, y)
	plt.plot(x, lr.predict(X), color='red', zorder=5)	
	plt.savefig('boostrapDemo.png')

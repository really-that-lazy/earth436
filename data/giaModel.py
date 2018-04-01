## giaModel.py #################################################################
## attempt to model the gia between sites using the data #######################
## in reformattedData.ods ######################################################
################################################################################
import pyexcel_ods

import sys

import random

import matplotlib.pyplot as plt
import numpy as np

from scipy import stats


from matplotlib.font_manager import FontProperties

import itertools

from rawhide import bootstrapper

from linearInterpolationModel import *
from giaUtils import *

fontP = FontProperties()
fontP.set_size('small')


def trendline(x, gradient, intercept):
	## return a y given an mx+b
	output = gradient*x + intercept		
	return output


## getLinearModel: listof(Num) listof(Num) -> listof(Num) Num Num Num listof(Num) listof(Num)

def getLinearModel(x_values, y_values, k=1.0, l=1.0):
	gradient, intercept, r_value, p_value, std_err = stats.linregress(x_values,y_values)	

	y_model = []
	yModelHigh = []
	yModelLow = []
	
	grad = k*gradient
	interc = l*intercept
	
	for x in x_values:
		y = trendline(x, grad, interc)
		yHigh = trendline(x, grad+(1.96*std_err), interc)
		yLow = trendline(x, grad-(1.96*std_err), interc)		
		y_model.append(y)
		yModelHigh.append(yHigh)
		yModelLow.append(yLow)	
	
	rSquare = r_value**2
	
	return y_model, grad, interc, std_err, yModelHigh, yModelLow, rSquare
	
	## yModelHigh and yModelLow are the y model built with a slope at the
	## extreme of the error bounds on the gradient



	

## ideas for future models:
## -same connect the dots idea, but with binned means every so many years
## -the extensions off the end use the means for both parameters in place of
## second largest/smallest age data point, as that part is royally fucking up
## the forecasts relative to the general trend










def plotGradientConfidenceIntervals(giaRegressionsByCombo, keys):
	def plotInterval(y, xstart, xstop, intervalLabel, colord, colords):
		"""Plot interval at y from xstart to xstop with given color."""   
  
		plt.hlines(y, xstart, xstop, colords, lw=7)
		plt.hlines(y, xstart, xstop, colord, lw=3, label=intervalLabel)
		##plt.vlines(xstart, y+0.3, y-0.3, colords, lw=2)
		##plt.vlines(xstop, y+0.3, y-0.3, colords, lw=2)


	y = 12
	for combo in keys:
		y -= 1
		combo1 = combo.split('-')[0]
		combo2 = combo.split('-')[1].split('_')[0]
		order = combo.split('-')[1].split('_')[1]
		
		if(order == 'forward'):
			d = combo1
			ds = combo2
		else:
			d = combo2
			ds = combo1
		if(order == 'forward'):
			plotInterval(y, abs(giaRegressionsByCombo[combo]['gradient'][0]), abs(giaRegressionsByCombo[combo]['gradient'][1]), "", mapSiteToColour(d), mapSiteToColour(ds))
		else:
			plotInterval(y, abs(giaRegressionsByCombo[combo]['gradient'][0]), abs(giaRegressionsByCombo[combo]['gradient'][1]), "",  mapSiteToColour(d), mapSiteToColour(ds))			
		
		est = abs(giaRegressionsByCombo[combo]['gradientEstimator'])
		
		plt.vlines(est, y+0.3, y-0.3, mapSiteToColour(d), lw=4)
	plt.xlabel('GIA')


	##plt.legend(loc=3, prop={'size': 11})
	plt.savefig('intervals.png')
	plt.close()




## "./reformattedData.ods"

def getDatasetsModelsAndObjects(filenameToLoad):
	lookupTable = pyexcel_ods.get_data(filenameToLoad)
	## open up the excel file to get the data as a dict of 2-lists
	locations = ['BATB', 'TAHB', 'GTB', 'ATB']
	## the first key for the lookupTable is the site location
	
	datasets = {}
	
	for loc in locations:
		datasets[loc] = [row for row in lookupTable[loc]]
		## under each key is a rectangular list with two columns to each row, 
		## the first one is elevation, the second one is age
		
	for d in datasets:
		print d, datasets[d], "\n\n\n"
	
	
	datasetObjects = {}
	
	datasetModels = {}
	
	for d in datasets:
		datasetObjects[d] = siteData(d, datasets[d])
		## build the dataset containers using the data retrieved for each site
		
		## note that the siteData object automatically filters the data received
		## to get rid of the first few non data lines and any empty spaces

	

	return datasets, datasetModels, datasetObjects



if(__name__ == "__main__"):

	datasets, datasetModels, datasetObjects = getDatasetsModelsAndObjects("./reformattedData.ods")
	
	allAgesSampled = [datasetObjects[d].getAgeValues() for d in datasets]
	allAgesSampled = [item for sublist in allAgesSampled for item in sublist]
	## flatten out the 2-list with some list comprehension
	print min(allAgesSampled), max(allAgesSampled)

	## create the raw plot of data points ######################################
	for ds in datasetObjects:
		print ds, datasetObjects[ds].data
		x = datasetObjects[ds].getAgeValues()
		y = datasetObjects[ds].getElevationValues()
		plt.plot(x, y, mapSiteToColour(ds) + 's', label=ds, markersize=4.0)
		
		datasetModels[ds] = siteModelConnectTheDots(datasetObjects[ds])
	
	plt.title("Plot of Elevation by Age\nRaw Data only")
	plt.ylabel('Elevation')
	plt.xlabel('Age')
	plt.legend(loc=2, prop={'size': 17})
	plt.savefig('./theDataRaw.png')
	plt.close()
	############################################################################
	




	## create the raw plot with the model included #############################
	for ds in datasetObjects:
		print ds, datasetObjects[ds].data
		x = datasetObjects[ds].getAgeValues()
		y = datasetObjects[ds].getElevationValues()
		plt.plot(x, y, mapSiteToColour(ds) + 's', label="%s, n=%i" % (ds, len(x)), markersize=4.0)

	for d in datasets:
		plt.plot(sorted(allAgesSampled), [datasetModels[d].getModelledElevation(age) for age in sorted(allAgesSampled)], mapSiteToColour(d), label=d+" (model)")

	plt.title("Plot of Elevation by Age\nRaw with Model")
	plt.ylabel('Elevation')
	plt.xlabel('Age')
	plt.legend(loc=2, prop={'size': 17})
	plt.savefig('./theData.png')
	plt.close()
	############################################################################


	############################################################################
	## create a list of the shortform name of all the sites, then use ##########
	## itertools to make a list of the sites!, all possible combinations of ####
	## sites ###################################################################
	
	
	## ie [A,B,C] -> [[A,B], [B,C], [C, A]]
	sites = [ds for ds in datasetObjects]
	siteCombinations = list(itertools.combinations(sites, 2))
	############################################################################
	
	for combo in siteCombinations:
		print combo
		for i in range(len(combo)):
			print i, combo[i]
	
	
	
	
	giaRegressions = {}
	
	giaRegressionsByCombo = {}
	
	giaRegressionKeys = []
	
	globalHistogramFloor = None
	## ayy lmao
	
	histogramFloorsList = []

	histogramFloorsByCombo = {}

	totalAges = []

	for combo in siteCombinations:
		
		histogramFloor = None
		ageFloor = None
		for site in combo:
			x = datasetObjects[site].getAgeValues()
			totalAges += x
			y = datasetObjects[site].getElevationValues()
			if(histogramFloor == None):
				histogramFloor = min(y)
			else:
				histogramFloor = min([min(y), histogramFloor])	

			if(ageFloor == None):
				ageFloor = min(x)
			else:
				ageFloor = min([min(x), ageFloor])

		def roundFloatDownToNearestTen(someFloat):
			someFloat /= 10
			someFloat = int(someFloat)
			someFloat *= 10
			someFloat -= 10
			return someFloat
					
		histogramFloor = roundFloatDownToNearestTen(histogramFloor)
		ageFloor = roundFloatDownToNearestTen(ageFloor)
		histogramFloorsByCombo[combo] = histogramFloor
			
		histogramFloorsList.append(ageFloor)
		print "histogramFloor for site combo", combo, ": ", histogramFloor
	
	globalHistogramFloor = min(histogramFloorsList)	
	
	print "global bin floor set at ", globalHistogramFloor
	
	globalBins=range(globalHistogramFloor, int(max(totalAges))+200, 200)
	
	print "global bins: ", globalBins
	
	for combo in siteCombinations:
		for i in globalBins:
			for site in combo:
				thisSiteDataset = datasetObjects[site]
				print "(",i, ",", i+200, "), site %.4s  count: " % site, thisSiteDataset.getThisSiteBinCount(i, 200)
	
	for combo in siteCombinations:

		
		
		for site in combo:
			
			
			
			print site 
			##datasetObjects[site].data
			x = datasetObjects[site].getAgeValues()
			y = datasetObjects[site].getElevationValues()

			plt.hist(x, bottom = histogramFloor, normed=False, bins=globalBins, alpha=0.4, color=mapSiteToColour(site))

			plt.plot(x, y, mapSiteToColour(site) + 's', label="%s, n=%i" % (site, len(x)), markersize=4.0)

			plt.plot(sorted(allAgesSampled), [datasetModels[site].getModelledElevation(age) for age in sorted(allAgesSampled)], mapSiteToColour(site), label=site+" (model)")



		plt.title("Data and Model for site Combination %s/%s" % (combo[0], combo[1]))	
		plt.ylabel('Elevation')
		plt.xlabel('Age')
		plt.legend(loc=2, prop={'size': 17})
		
		
		axes1 = plt.gca()
		print plt.gca()
		yScaleRange =  max(axes1.get_ylim()) - min(axes1.get_ylim())


		axes2 = plt.twinx()
		axes2.set_ylabel('Count')		
		axes2.axis([None,None,0,yScaleRange])
		
		
		plt.savefig('%s-%s_DataAndModel.png' % (combo[0], combo[1]))
		
		plt.close()
	
	
		for order in ['forward', 'reverse']:
			if(order == 'forward'):
				d = combo[0]
				## d is the site we are using as our direct comparison
				
				## ie MUST have a measured data point at this age
				ds = combo[1]
				## ds is what we are comparing against, so it can just be a
				## modelled point
			else:
				d = combo[1]
				ds = combo[0]
			
	
			
			allowableAgeValues = []
			
			for age in sorted(allAgesSampled):
				if(datasetModels[d].ageValueInRawData(age) and datasetModels[ds].ageValueIsInRangeCoveredByModel(age) and datasetModels[ds].ageComparisonValidForThisBin(datasetModels[d], globalBins, age) ):
					allowableAgeValues.append(age)
				else:
					continue
					## the case where we have an overlap of the models, but
					## either A: no datapoint is actually present for either
					## dataset at this age, so comparisons are not honouring
					## the raw data, or
					## B: we have a datapoint on the set to compare against
					## but not the one we are comparing
					##else:
						##continue
						## if the datapoint in question is outside the bounds
						## covered by these two datasets, they cant be considered
			
			##allowableAgeValues = [age for age in sorted(allAgesSampled) if( (datasetModels[ds].ageValueIsInRangeCoveredByModel(age))and(datasetModels[d].ageValueIsInRangeCoveredByModel(age))and( (datasetModels[d].ageValueInRawData(age))or(datasetModels[ds].ageValueInRawData(age)) )   )]
			
			elevationDiffs = [(datasetModels[d].getModelledElevation(age) - datasetModels[ds].getModelledElevation(age)) for age in allowableAgeValues]
			
			thisComparisonGia = " measured data from %s relative to %s model" % (d, ds)
	
			bootstrapper.plotBootstrapsOnDataPlot(plt, allowableAgeValues, elevationDiffs, mapSiteToColour(ds), mapSiteToColour(d));		

			
			plt.plot(allowableAgeValues, elevationDiffs, mapSiteToColour(d)+'+', label=thisComparisonGia, markersize=4.0)
			
			
			linRegressYValues, gradient, intercept, gradientError, 	yModelHigh, yModelLow, rSquare = getLinearModel(allowableAgeValues, elevationDiffs)
			
			if(d != ds):
				giaRegressions[thisComparisonGia] = {"N": len(allowableAgeValues), "gradientEstimator": gradient, "gradientError": gradientError, "gradient": [gradient+(1.96*gradientError), gradient-(1.96*gradientError)], "intercept": intercept, }
				giaRegressionsByCombo["%s-%s_%s" % (combo[0], combo[1], order)] = {"N": len(allowableAgeValues), "gradientEstimator": gradient, "gradientError": gradientError, "gradient": [gradient+(1.96*gradientError), gradient-(1.96*gradientError)], "intercept": intercept, }			
				giaRegressionKeys.append("%s-%s_%s" % (combo[0], combo[1], order))
			

			##plt.plot(allowableAgeValues, linRegressYValues, mapSiteToColour(d)+'', label=" measured data from %s relative to %s model linearRegression" % (d, ds), markersize=4.0)
			##plt.plot(allowableAgeValues, yModelHigh, mapSiteToColour(ds)+'', markersize=4.0)
			##plt.plot(allowableAgeValues, yModelLow, mapSiteToColour(ds)+'', markersize=4.0)
			
	
			plt.suptitle("Plot of Elevation Diff for %s relative to %s model by Age\n(%s order for %s-%s), n=%i\ny=mx+b, m = %.4f SE(%.4f), b = %.4f, r^2 = %.3f" % (d, ds, order, combo[0], combo[1], len(allowableAgeValues), gradient, gradientError, intercept, rSquare), fontsize=10)
			plt.ylabel('Elevation')
			plt.xlabel('Age')
			if(d == "ATB"):
				plt.legend(loc=2, prop={'size': 14})

			else:	
				plt.legend( loc=3, prop={'size': 14})
			##plt.savefig('./theGIA_%s_relative_to_%s.png' % (d, ds))
			## ^ this was creating a ton of clutter
			plt.savefig('./gias/theGIA_%s_relative_to_%s.png' % (d, ds))
			plt.close()
	
	
	
	ageMatches = []
	for d in datasets:
		for dv in datasetObjects[d].getAgeValues():
			for od in datasets:
				if(od != d):
					if((dv in datasetObjects[od].getAgeValues())and dv not in ageMatches):
						ageMatches.append(dv)
	print ageMatches
	
	sortedKeys = sorted([regress for regress in giaRegressions])
	
	for regress in sortedKeys:
		print regress			
	
	print "\n\n"
	
	for regress in sortedKeys:
		print regress, ": ", giaRegressions[regress], "\n"
	
	print "\n\n\n"		

	

	for regress in sortedKeys:
		print regress, ",", giaRegressions[regress]['gradient'][0], ",", giaRegressions[regress]['gradient'][1]	
	
	plotGradientConfidenceIntervals(giaRegressionsByCombo, giaRegressionKeys)

	for site in sites:
		plt.plot([1], [1], mapSiteToColour(site)+'s', label=site, markersize=20)
		plt.plot([1], [1], mapSiteToColour(site), label=site+" model", markersize=20)
	plt.axis('off')
	plt.legend(loc=3, prop={'size': 29})
	plt.savefig("legendary.png")
	plt.close()

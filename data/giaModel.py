## giaModel.py #################################################################
## attempt to model the gia between sites using the data #######################
## in reformattedData.ods ######################################################
################################################################################
import pyexcel_ods

##import sys


import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from matplotlib.font_manager import FontProperties

import itertools
## used to generate the number of links between sites
import random
from random import sample
## used in the random choice feature of the zoomed site comparisons 

from rawhide import bootstrapper
## get the custom bootstrap plotting function from this projects code

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



	












def plotGradientConfidenceIntervals(giaRegressionsByCombo, keys, giaRegressionDescriptions, outputPathDict):
	def plotInterval(ax, y, xstart, xstop, intervalLabel, colord, colords):
		"""Plot interval at y from xstart to xstop with given color."""   
  
		ax.hlines(y, xstart, xstop, colords, lw=7)
		ax.hlines(y, xstart, xstop, colord, lw=3, label=intervalLabel)
		## plots the interval in the colours of both sites


	outputPath = convertListToRelativePath([outputPathDict[setting] for setting in getCurrentSettingOptions()])
	
	

	y = 0
	## used in spacing out the intervals for each site vertically through the
	## graph
	
	fig,ax = plt.subplots(1)
	
	for combo in keys:
		y += 1
		combo1 = combo.split('-')[0]
		combo2 = combo.split('-')[1].split('_')[0]
		order = combo.split('-')[1].split('_')[1]
		
		if(order == 'forward'):
			direct = combo1
			modelled = combo2
		else:
			direct = combo2
			modelled = combo1
			
		est = giaRegressionsByCombo[combo]['gradientEstimator']	
			
		ciStart = giaRegressionsByCombo[combo]['gradient'][0]
		ciEnd = giaRegressionsByCombo[combo]['gradient'][1]
		
		if(est < 0):
			est = -est
			ciStart = -ciStart
			ciEnd = -ciEnd
			
		if(order == 'forward'):
			plotInterval(ax, y, ciStart, ciEnd, "", mapSiteToColour(direct), mapSiteToColour(modelled))
		else:
			plotInterval(ax, y, ciStart, ciEnd, "",  mapSiteToColour(direct), mapSiteToColour(modelled))			
		
		
		ax.vlines(est, y+0.3, y-0.3, mapSiteToColour(direct), lw=4)
	ax.set_xlabel('GIA (m/year)')
	
	ax.set_xlim([0,0.009])

	plt.yticks(list(np.arange(1, len(keys)+1, 1.0)), [giaRegressionDescriptions[key] for key in keys], rotation=0)

	fileNameIdentifier = "_".join([outputPathDict[setting] for setting in getCurrentSettingOptions()])

	plt.title("95p Confidence intervals on GIA\nfilters: %s" % fileNameIdentifier)

	for item in ax.get_yticklabels():
		item.set_fontsize(8)
	
	outputFilePath = filePathOnRelativePath(outputPath+"gias/", fileName='intervals', ext="png")
	print "Saving gia intervals plot at '%s'" % outputFilePath
	verifyPath(outputPath+"gias/")
	

	
	plt.savefig(outputFilePath,bbox_inches='tight')





	outputFilePath = filePathOnRelativePath(outputPath+"gias/", fileName='%s_intervals' % fileNameIdentifier, ext="png")

	print "Saving gia intervals plot at '%s'" % outputFilePath

	plt.savefig(outputFilePath,bbox_inches='tight')

	plt.close()




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
	for site in datasetObjects:
		print site, datasetObjects[site].data
		x = datasetObjects[site].getAgeValues()
		y = datasetObjects[site].getElevationValues()
		
		n = len(datasetObjects[site].getAgeValues())
		
		plt.plot(x, y, mapSiteToColour(site) + 's', label=site+" n=%i" % n, markersize=4.0)
		
		datasetModels[site] = siteModelConnectTheDots(datasetObjects[site])
	
	##plt.title("Plot of Elevation by Age\nRaw Data only")
	plt.ylabel('Elevation (m IGLD1985)')
	plt.xlabel('Age Before Present (years)')
	plt.legend(loc=2, prop={'size': 17})
	plt.savefig('./theDataRaw.png')
	plt.close()
	############################################################################
	



	############################################################################
	## create the raw plot with the model included #############################
	for ds in datasetObjects:
		print ds, datasetObjects[ds].data
		x = datasetObjects[ds].getAgeValues()
		y = datasetObjects[ds].getElevationValues()
		plt.plot(x, y, mapSiteToColour(ds) + 's', label="%s, n=%i" % (ds, len(x)), markersize=4.0)

	for d in datasets:
		plt.plot([age for age in sorted(allAgesSampled)  if datasetModels[d].ageValueIsInRangeCoveredByModel(age)], [datasetModels[d].getModelledElevation(age) for age in sorted(allAgesSampled) if datasetModels[d].ageValueIsInRangeCoveredByModel(age)], mapSiteToColour(d), label=d+" (model)")

	##plt.title("Plot of Elevation by Age\nRaw Data with Model")
	plt.ylabel('Elevation (m IGLD1985)')
	plt.xlabel('Age Before Present (years)')
	plt.legend(loc=2, prop={'size': 17})
	plt.savefig('./theData.png')
	plt.close()
	############################################################################
	
	############################################################################
	## create the raw plot with the model included, zooming in on the 2000-2300#
	## ybp window, y axis limited to 183-187 m #################################
	
	
	zoomXRange = (2000, 2400)
	zoomYRange = (182, 190)
	
	for site in datasetObjects:
		print site, datasetObjects[site].data
		x = datasetObjects[site].getAgeValues()
		y = datasetObjects[site].getElevationValues()
		plt.plot(x, y, mapSiteToColour(site) + 's', label="%s" % (site), markersize=4.0)

	for site in datasets:
		plt.plot(sorted(allAgesSampled), [datasetModels[site].getModelledElevation(age) for age in sorted(allAgesSampled)], mapSiteToColour(site), label=site+" (model)")
		## plot the dataset models as straight lines
		
	siteCodeOptions = [site for site in datasetObjects]
	
	exampleSites = sample(siteCodeOptions, 2)
	print exampleSites
	
	direct = exampleSites[0]
	modelled = exampleSites[1]
	
	agesToConsider = [age for age in sorted(allAgesSampled) if ( ((age >= min(zoomXRange))and(age <= max(zoomXRange))) and (datasetModels[direct].ageValueInRawData(age) and datasetModels[modelled].ageValueIsInRangeCoveredByModel(age)) )]
	
	for age in agesToConsider:
		print age
	
	demoComparisonPoint = random.choice(agesToConsider)
	print "-> ", demoComparisonPoint	
	
	
	directElevation = datasetObjects[direct].getElevationByGivenAge(demoComparisonPoint)
	modelledElevation = datasetModels[modelled].getModelledElevation(demoComparisonPoint)	
	
	print "Direct [%s]:" % direct, directElevation
	print "Modelled [%s]:" % modelled, modelledElevation
	
	##exit()
	##ax = plt.axes()
	##ax.arrow(demoComparisonPoint, directElevation, 0, modelledElevation-directElevation, head_width=5.5, head_length=10.1, fc=mapSiteToColour(direct), ec="y")
	plt.plot([demoComparisonPoint, demoComparisonPoint], [directElevation, modelledElevation], "%s" % mapSiteToColour(direct), linewidth=3.0)
	plt.plot([demoComparisonPoint, demoComparisonPoint], [directElevation, modelledElevation], "%s--" % mapSiteToColour(modelled), linewidth=2.0)
	##'--'
	plt.title("Plot of Elevation by Age\nRaw Data with Model")
	plt.ylabel('Elevation (m)')
	plt.xlabel('Age Before Present (years)')
	plt.axis((zoomXRange[0], zoomXRange[1],zoomYRange[0],zoomYRange[1]))
	plt.legend(loc=2, prop={'size': 7})
	plt.savefig('./theDataZoomed.png')
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
	
	
	
	

	
	
	
	globalHistogramFloor = None
	## ayy lmao
	
	histogramFloorsList = []

	histogramFloorsByCombo = {}

	totalAges = []



	############################################################################
	## loop through the site combinations and use the data to decide on a ######
	## floor for the age bins and the bounds on the plot axes ##################
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
	############################################################################

	
	globalHistogramFloor = min(histogramFloorsList)	
	
	print "global bin floor set at ", globalHistogramFloor
	
	globalBins=range(globalHistogramFloor, int(max(totalAges))+200, 200)
	## build a list of bin endpoints starting at the floor value and ending at
	## one bin width above the last age value of any of the dataset
	
	## example of how this works if you run
	## range(450, 4857+200, 200)	
	print "global bins: ", globalBins
	
	############################################################################
	## for debug output print out bin counts for each dataset ##################
	for i in globalBins:
		print "bin (",i, ",", i+200, "):"
		print "-"*80
		for combo in siteCombinations:
		
			for site in combo:
				thisSiteDataset = datasetObjects[site]
				
				siteName = '{:4s}'.format(site)
				
				print "site %.4s  count: %i" % (siteName, thisSiteDataset.getThisSiteBinCount(i, 200))
			print "-"*80
		print "\n\n"
	############################################################################
	
	
	giaRegressions = {}
	giaRegressionComboMappingsByConditions = {}	
	
	giaRegressionKeys = []
	giaRegressionDescriptions = {}
	giaKeysByDescriptions = {}
	
	for combo in siteCombinations:
		for order in ['forward', 'reverse']:
			if(order == 'forward'):
				direct = combo[0]
				## d is the site we are using as our direct comparison
				
				## ie MUST have a measured data point at this age
				modelled = combo[1]
				## ds is what we are comparing against, so it can just be a
				## modelled point
			else:
				direct = combo[1]
				modelled = combo[0]
			thisRegressionKey = "%s-%s_%s" % (combo[0], combo[1], order)
				
			giaRegressionKeys.append(thisRegressionKey)
			thisComparisonGiaDescription = " measured data from %s relative to %s model" % (direct, modelled)
			giaRegressionDescriptions[thisRegressionKey] = thisComparisonGiaDescription
			giaKeysByDescriptions[thisComparisonGiaDescription] = thisRegressionKey
	
	sortedKeys =  sorted(giaRegressionKeys)
	print "sortedKeys: ", sortedKeys	

	############################################################################
	## plot the raw data plots with counts for each bin ########################
		
	for combo in siteCombinations:
		print "\nPlotting raw data for site combo: "
		for site in combo:
			print site		
			
		for site in combo: 
			x = datasetObjects[site].getAgeValues()
			y = datasetObjects[site].getElevationValues()


			plt.plot(x, y, mapSiteToColour(site) + 's', label="%s, n=%i" % (site, len(x)), markersize=4.0)
			## plot the raw data for each site
				
			plt.plot(sorted(allAgesSampled), [datasetModels[site].getModelledElevation(age) for age in sorted(allAgesSampled)], mapSiteToColour(site), label=site+" (model)")
			## plot the linear interpolation model for each site
				
			plt.hist(x, bottom = histogramFloor, normed=False, bins=globalBins, alpha=0.4, color=mapSiteToColour(site))
			## plot the histogram of data set counts on the plot alongside the 
			## data itself
				
			## histogram floor was chosen here as a nice looking spot to put the
			## count histogram so it doesnt overlap the main data


		plt.title("Data and Model for site Combination %s/%s" % (combo[0], combo[1]))	
		plt.ylabel('Elevation')
		plt.xlabel('Age')
		plt.legend(loc=2, prop={'size': 17})
			
			
		axes1 = plt.gca()
		yScaleRange =  max(axes1.get_ylim()) - min(axes1.get_ylim())


		axes2 = plt.twinx()
		axes2.set_ylabel('Count')		
		axes2.axis([None,None,0,yScaleRange])
		## set the left axis to be elevation relative to datum
			
		## and the right axis to be count of each dataset in each bin
			
		outputFilePath = filePathOnRelativePath("./", fileName='%s-%s_DataAndModel' % (combo[0], combo[1]), ext="png")
		print "Saving rawData plot at '%s'" % outputFilePath
		verifyPath("./")
		## umm ok
		plt.savefig(outputFilePath)
		plt.close()
		## save the raw data combo graph
	############################################################################	
		
	
	############################################################################
	## plot the gia graphs and store the raw regression numbers used to create #
	## them ####################################################################
	for conditions in [{"valueDifference": "withinThirtyPercent","valueCounts": "bothNonZero"},\
						{"valueDifference": "withinFiftyPercent","valueCounts": "bothNonZero"}, \
						##{"valueDifference": "withinTwentyPercent","valueCounts": "bothNonZero"}, \
						{"valueDifference": "withinSeventyFivePercent","valueCounts": "bothNonZero"},  \
						##{"valueDifference": "withinSixtyPercent","valueCounts": "bothNonZero"},  \
						{"valueCounts": "bothNonZero"}]:
		outputPathDict = populateConditionsDict(conditions)

		outputPath = convertListToRelativePath([outputPathDict[setting] for setting in getCurrentSettingOptions()])
		conditionIdString = "_".join([outputPathDict[setting] for setting in getCurrentSettingOptions()])
		
		giaRegressionsByCombo = {}
		
		for combo in siteCombinations:


			print "\nPlotting gia for site combo: "
			for site in combo:
				print site		
			
		
		
			## now the gia calculations
			for order in ['forward', 'reverse']:
				## each comparison has a forward A to B, and reverse B to A,
				## comparison, the CIs on the absolute value of slope for this must
				## be statistically similar for the comparison to work 
				
				if(order == 'forward'):
					direct = combo[0]
					## d is the site we are using as our direct comparison
					
					## ie MUST have a measured data point at this age
					modelled = combo[1]
					## ds is what we are comparing against, so it can just be a
					## modelled point
				else:
					direct = combo[1]
					modelled = combo[0]
				
		
				
				allowableAgeValues = []
				## for each comparison, there are only a small number of data values
				## from the initial dataset that can be used for valid comparison
				
				## each datapoint used for a gia comparison must be:
				## -from the direct dataset
				## -in the range covered by the modelled dataset (meaning that if 
				## the direct comparison dataset has a datapoint available, but the
				## modelled one has just been hanging off the end in a straight line
				## from the last known datapoint, it cant be considered valid
				## -given that theres a bin from startAge to startAge+binWidth that
				## the datapoints age is in, that bin needs to hit some criteria for
				## the number of datapoints in the bin from both  
				
				for age in sorted(allAgesSampled):
					if(datasetModels[direct].ageValueInRawData(age) and datasetModels[modelled].ageValueIsInRangeCoveredByModel(age) and datasetModels[modelled].ageComparisonValidForThisBin(datasetModels[direct], globalBins, age, conditions) ):
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
				

				elevationDiffs = [(datasetModels[direct].getModelledElevation(age) - datasetModels[modelled].getModelledElevation(age)) for age in allowableAgeValues]
				
		
				bootstrapper.plotBootstrapsOnDataPlot(plt, allowableAgeValues, elevationDiffs, mapSiteToColour(modelled), mapSiteToColour(direct));		

				thisComparisonGiaDescription = " measured data from %s relative to %s model" % (direct, modelled)
				plt.plot(allowableAgeValues, elevationDiffs, mapSiteToColour(direct)+'+', label=thisComparisonGiaDescription, markersize=4.0)
				
				
				linRegressYValues, gradient, intercept, gradientError, 	yModelHigh, yModelLow, rSquare = getLinearModel(allowableAgeValues, elevationDiffs)
				
				if(direct != modelled):
					giaRegressionKey = "%s-%s_%s" % (combo[0], combo[1], order)
					
					giaRegressionsByCombo[giaRegressionKey] = {"N": len(allowableAgeValues), "gradientEstimator": gradient, "gradientError": gradientError, "gradient": [gradient+(1.96*gradientError), gradient-(1.96*gradientError)], "intercept": intercept, }			
		
				plt.suptitle("Plot of Elevation Diff for %s relative to %s model by Age\n(%s order for %s-%s), n=%i\ny=mx+b, m = %.4f SE(%.4f), b = %.4f, r^2 = %.3f" % (direct, modelled, order, combo[0], combo[1], len(allowableAgeValues), gradient, gradientError, intercept, rSquare), fontsize=10)
				plt.ylabel('Elevation')
				plt.xlabel('Age')
				if(direct == "ATB"):
					plt.legend(loc=2, prop={'size': 14})

				else:	
					plt.legend( loc=3, prop={'size': 14})
				##plt.savefig('./theGIA_%s_relative_to_%s.png' % (d, ds))
				## ^ this was creating a ton of clutter


				outputFilePath = filePathOnRelativePath(outputPath+"gias/", fileName='theGIA_%s_relative_to_%s' % (direct, modelled), ext="png")
				print "Saving gia plot at '%s'" % outputFilePath
				verifyPath(outputPath+"gias/")
				plt.savefig(outputFilePath)
				plt.close()
		giaRegressionComboMappingsByConditions[conditionIdString] = giaRegressionsByCombo
		plotGradientConfidenceIntervals(giaRegressionsByCombo, giaRegressionKeys, giaRegressionDescriptions, outputPathDict)
	############################################################################
	print "Finished gia plots"


	############################################################################
	## Check for any exact age matches in the datasets provided ################
	## Spoiler: there arent any ################################################
	ageMatches = []
	for d in datasets:
		for dv in datasetObjects[d].getAgeValues():
			for od in datasets:
				if(od != d):
					if((dv in datasetObjects[od].getAgeValues())and dv not in ageMatches):
						ageMatches.append(dv)
	print "Exact age matches between datasets: ", ageMatches
	############################################################################
	



	############################################################################
	## now that values have been generated for GIA for each site comparison, ###
	## convert them to intervals for each site combination and save the result #
	## to file #################################################################
	for idString in giaRegressionComboMappingsByConditions:
		print "\n\n%s:	" % idString
		
		giaRegressionsByCombo = giaRegressionComboMappingsByConditions[idString]
		siteCombos = ["ATB-BATB","GTB-ATB","GTB-BATB","GTB-TAHB","TAHB-ATB","TAHB-BATB"]

		with open("%s_intervals.csv" % idString, "wb") as csv_file:
			writer = csv.writer(csv_file, delimiter=',')
			
			writer.writerow(["id", "name", "startValue", "endValue"])
			for regress in sortedKeys:
				description = giaRegressionDescriptions[regress]
				ciStart = giaRegressionsByCombo[regress]['gradient'][0]
				ciEnd = giaRegressionsByCombo[regress]['gradient'][1]
				

				writer.writerow([regress, description, ciStart, ciEnd])
		
			for combo in siteCombos:
				for order in ["forward", "reverse"]:
					regress = "%s_%s" % (combo, order)
					print giaRegressionDescriptions[regress], ",", giaRegressionsByCombo[regress]['gradient'][0], ",", giaRegressionsByCombo[regress]['gradient'][1]			
				
					est = giaRegressionsByCombo[regress]['gradientEstimator']	
					ciStart = 100*100*giaRegressionsByCombo[regress]['gradient'][0]
					ciEnd = 100*100*giaRegressionsByCombo[regress]['gradient'][1]
				
					if(est < 0):
						ciStart = -ciStart
						ciEnd = -ciEnd
					
					if(order == "forward"):
						forwardInterval = {"start":min(ciStart, ciEnd), "end": max(ciStart, ciEnd)}
					elif(order == "reverse"):
						reverseInterval = {"start":min(ciStart, ciEnd), "end": max(ciStart, ciEnd)}
				mergedInterval = mergeConfidenceIntervals(forwardInterval, reverseInterval)
				
				print combo, ": ", mergedInterval
				if(mergedInterval == "No overlap"):
					writer.writerow([combo, "%s_merged" % combo, mergedInterval, ""])
				else:
					writer.writerow([combo, "%s_merged" % combo, "%.3f" % mergedInterval[0], "%.3f" % mergedInterval[1]])

	

	############################################################################
	## plot a legend showing the colour coding system for the sites ############
	## this sounded like a decent idea earlier, but it eventually proved #######
	## not to be needed ########################################################
	for site in sites:
		plt.plot([1], [1], mapSiteToColour(site)+'s', label=site, markersize=20)
		plt.plot([1], [1], mapSiteToColour(site), label=site+" model", markersize=20)
	plt.axis('off')
	plt.legend(loc=3, prop={'size': 29})
	plt.savefig("legendary.png")
	plt.close()
	############################################################################

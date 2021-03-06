## linearInterpolationModel.py #################################################
## model for elevation v. time values made by connecting point to point, ie ####
## "connect the dots" formally known as linear interpolation ###################
################################################################################
from dataModel import *

import numpy as np


def percentageDifference(someValue, anotherValue):
	average = float(someValue + anotherValue)/2.0
	diff = abs(someValue - anotherValue)
	return float(diff/average)
	
	
def conditionMet(thisBinCount, otherBinCount, condition):
	if(condition == "any"):
		return True
	elif(condition == "bothNonZero"):
		if((thisBinCount > 0) and (otherBinCount > 0)):
			return True
		return False
	elif(condition == "withinTwentyPercent"):
		if(percentageDifference(thisBinCount, otherBinCount) <= 0.20):
			return True
		return False
	elif(condition == "withinThirtyPercent"):
		if(percentageDifference(thisBinCount, otherBinCount) <= 0.30):
			return True
		return False
	elif(condition == "withinFiftyPercent"):
		if(percentageDifference(thisBinCount, otherBinCount) <= 0.50):
			return True
		return False
	elif(condition == "withinSixtyPercent"):
		if(percentageDifference(thisBinCount, otherBinCount) <= 0.60):
			return True
		return False	
	elif(condition == "withinSeventyFivePercent"):
		if(percentageDifference(thisBinCount, otherBinCount) <= 0.75):
			return True
		return False		

	
class siteModelConnectTheDots(siteModel):	
	
	## siteModelConnectTheDots: siteData -> siteModelConnectTheDots
	
	def __init__(self, availableData):
		self.siteName = availableData.getSiteName()
		self.rawDataObject = availableData
	
	def ageValueIsInRangeCoveredByModel(self, someAge):
		if(someAge in self.rawDataObject.getAgeValues()):
			return True
		else:
			maxAgeCovered = max(self.rawDataObject.getAgeValues())
			minAgeCovered = min(self.rawDataObject.getAgeValues())
			if((someAge <= maxAgeCovered)and(someAge >= minAgeCovered)):
				return True
			else:
				return False
		
	def ageValueInRawData(self, someAge):
		if(someAge in self.rawDataObject.getAgeValues()):
			return True	
		return False

	def ageComparisonValidForThisBin(self, otherModelToCompareAgainst, globalBins, ageValue, conditions):
		binStart, binWidth = getAgeBinByAgeValue(ageValue, globalBins)
			
		thisModelsBinCount = self.rawDataObject.getThisSiteBinCount(binStart, binWidth)
		otherModelsBinCount = otherModelToCompareAgainst.rawDataObject.getThisSiteBinCount(binStart, binWidth)
		
		for condition in conditions:
			if(not conditionMet(thisModelsBinCount, otherModelsBinCount, conditions[condition])):
				return False
		return True
		## the loop successfully met every condition, so we're good to go
		
		
	def getModelledElevation(self, someAge):
		if(someAge in self.rawDataObject.getAgeValues()):
			return self.rawDataObject.getElevationByGivenAge(someAge)
		else:
			## dont have a datapoint available at that age value, so we need to
			## interpolate linearly between them to get it
			ageValues = np.array(self.rawDataObject.getAgeValues())

			if(ageValues[ageValues < someAge].size == 0):
				## case where our value to interpolate is off the bottom end
				## of the dataset, so we extrapolate from the last two values
				## min, and 2dmin
				
				minValue = min(ageValues)
				restOfValues = np.array([val for val in ageValues if val != minValue])
				## maybe npifying the array will make the min/max calls faster
				## idk
				secondMinValue = min(restOfValues)
				ageDelta = someAge - secondMinValue
				## distance from second smallest to the point we want to
				## interpolate
				
				
				
				secondMinAgeElevation = self.rawDataObject.getElevationByGivenAge(secondMinValue)
				minAgeElevation = self.rawDataObject.getElevationByGivenAge(minValue)
				
				outputElevationGuess = secondMinAgeElevation + ( (ageDelta/(abs(secondMinValue-minValue)))*(secondMinAgeElevation - minAgeElevation) )	
			
			elif(ageValues[ageValues > someAge].size == 0):
				## case where our value to interpolate is off the top end
				maxValue = max(ageValues)
				restOfValues = np.array([val for val in ageValues if val != maxValue])
				## npifying this array could make the min/max calls faster
				
				secondMaxValue = max(restOfValues)
				
				ageDelta = someAge - secondMaxValue
				
				secondMaxAgeElevation = self.rawDataObject.getElevationByGivenAge(secondMaxValue)
				maxAgeElevation = self.rawDataObject.getElevationByGivenAge(maxValue)				

				outputElevationGuess = secondMaxAgeElevation + ( (ageDelta/(abs(maxValue - secondMaxValue)))*(maxAgeElevation - secondMaxAgeElevation) )				

			else:
				closestAgeBelow = ageValues[ageValues < someAge].max()
				closestAgeAbove = ageValues[ageValues > someAge].min()			

				
			
			
				ageDelta = closestAgeAbove - closestAgeBelow
			
				elevBelow = self.rawDataObject.getElevationByGivenAge(closestAgeBelow)
				elevAbove = self.rawDataObject.getElevationByGivenAge(closestAgeAbove)
				
				elevDelta = elevAbove - elevBelow
			
				outputElevationGuess = elevBelow + elevDelta*((someAge - closestAgeBelow)/(ageDelta))
			
			return outputElevationGuess
		

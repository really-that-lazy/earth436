## dataModel.py ################################################################
## Base class for building data models (ie interpretations of what the age #####
## is for the entire age range that the data spans) ############################
################################################################################
from rawData import *


## ideas for future models:
## -same connect the dots idea, but with binned means every so many years
## -use mean water level of other sites in areas where a site has a gap to 
## adjust for global water lows caused by climate, etc.


class siteModel(object):
	## parent to all models that take a set of siteData and attempt to build a
	## model of elevations for all of the possible age values in between points
	## where the elevation is directly sampled for that exact time.
	
	
	
	def __init__(self):
		pass
	
	## getModelledElevation: Num -> Num




def inRange(value, high, low):
	if((value <= high) and (value >= low)):
		return True
	return False


## getAgeBinByAgeValue: float, listof(float) -> float, float

## ageValue is a float
## ageBins is some list of bin endpoints

def getAgeBinByAgeValue(ageValue, ageBins):
	baseAge = ageBins[0]
	ageBinsDelta = ageBins[1] - baseAge
	## should be consistent throughout the ageBins list
	
	## now find the nearest start point for a bin to the ageValue
	
	binStartValue = baseAge
	
	if(ageValue < baseAge):
		while(True):
			if(inRange(ageValue, binStartValue+ageBinsDelta, binStartValue)):
				return binStartValue, ageBinsDelta			
			binStartValue -= ageBinsDelta		
	else:
		while(True):
			if(inRange(ageValue, binStartValue+ageBinsDelta, binStartValue)):
				return binStartValue, ageBinsDelta			
			binStartValue += ageBinsDelta


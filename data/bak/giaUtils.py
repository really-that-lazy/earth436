## giaUtils.py #################################################################
## utility functions for gia thesis ############################################b
################################################################################
import os


def mapSiteToColour(siteLoc):
	siteMappings ={'BATB': 'b', 'TAHB': 'g', 'GTB': 'm', 'ATB': 'r'}
	return siteMappings[siteLoc]

def convertListToRelativePath(someListOfStrings):
	output = "./"
	for someStr in someListOfStrings:
		output += "%s/" % someStr
	return output

def filePathOnRelativePath(somePath, fileName, ext=None):
	if(ext == None):
		return "%s%s" % (somePath, fileName)
	else:
		return "%s%s.%s" % (somePath, fileName, ext)

def verifyPath(somePath):
	if(os.path.exists(somePath)):
		if(os.path.isdir(somePath)):
			return True
		else:
			print "Path '%s' exists, but is not a directory"
			return False
	else:
		print "Path did not exist, attempting to create it..."
		os.makedirs(somePath)
		return os.path.exists(somePath)

## ie {"valueDifference": "withinTwentyPercent","valueCounts": "bothNonZero"}
		
def populateConditionsDict(inputDict):
	if("valueDifference" not in inputDict):
		inputDict["valueDifference"] = "any"
	if("valueCounts" not in inputDict):
		inputDict["valueCounts"] = "any"
	return inputDict
	## technically it mutates the dict, but this is ok too

def getCurrentSettingOptions():
	return [ "valueCounts", "valueDifference"]
		
		
if(__name__ == "__main__"):
	print convertListToRelativePath(["withinTwentyPercent", "baseFixedAt450", "gias"])

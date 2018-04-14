## rawData.py ##################################################################
## object containing the raw data object used to wrap what gets pulled out #####
## of the spreadsheet file #####################################################
################################################################################



class siteData(object):
	
	## siteData: Str Listof(Any) -> siteData
	
	def __init__(self, siteName, rawData):
		self.dataHeader = []
		self.data = [row for row in rawData if len(row) == 2]
		## filter to only those rows that contain exactly two elements in order
		## to get rid of the unimportant details at the top of each sheet
		self.siteName = siteName
		
		
	def getAgeValues(self):
		## get all available (measured) data points for age
		return [row[1] for row in self.data]
		

	def getElevationValues(self):
		return [row[0] for row in self.data]
		
	def getElevationByGivenAge(self, someAge):	
		for row in self.data:
			if(row[1] == someAge):
				return row[0]

	def getThisSiteBinCount(self, binStart, binWidth):
		
		def withinside(someValue, binStart, binWidth):
			delta = someValue - binStart
			if((delta >= 0)and(delta <= binWidth)):
				return True
			else:
				return False
		
		return len([row[1] for row in self.data if withinside(row[1], binStart, binWidth)])	
		
	def getSiteName(self):
		return self.siteName

#import numpy as np
import csvData
import csvFree
from myUtils import getStarWithMinDist

lamostFile = '/Volumes/obiwan/azuri/data/lamost/dr7_v0_q1.csv'#'/Volumes/obiwan/azuri/data/lamost/dr7_med_v0_q1.csv'

lamostCSV = csvFree.readCSVFile(lamostFile, '|')

pa30_ra = 13.29667
pa30_dec = 67.50067

closestStars = getStarWithMinDist(lamostCSV, pa30_ra, pa30_dec)
print('closestStars = ',closestStars)

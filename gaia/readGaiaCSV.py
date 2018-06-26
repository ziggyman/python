import os
import gzip
import numpy as np

dataDir = '/Volumes/external/azuri/data/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/csv/'
fileNameRoot = 'GaiaSource_%03d-%03d-%03d.csv.gz'#% (release, tract, patch)

filename = os.path.join(dataDir, fileNameRoot % (0, 0, 0))
with gzip.open(filename, 'rb') as f:
    file_content = f.read()
    file_lines = file_content.splitlines()
    header = file_lines[0].split(',')
print 'header = ',header

fileLines = []
dataList = []
previous = False
prePrevious = False
isThere = False
for iRelease in range(3):
    for iTract in np.arange(0, 1000):
        for iPatch in np.arange(0, 1000):
            prePrevious = previous
            previous = isThere
            filename = os.path.join(dataDir, fileNameRoot % (iRelease, iTract, iPatch))
            if os.path.exists(filename):
                isThere = True
                if prePrevious and not previous:
                    print 'previous file is missing'
                with gzip.open(filename, 'rb') as f:
                    file_content = f.read()
                    file_lines = file_content.splitlines()
                    header = file_lines[0].split(',')
                    del file_lines[0]
    #                fileLines.append(file_lines)
                    for line in file_lines:
                        line_split = line.split(',')
#                        print 'line_split = ',line_split
                        data = {}
                        for iCol in range(len(header)):
                            data[header[iCol]] = line_split[iCol]
#                            print 'data[',header[iCol],'] = ',data[header[iCol]]
    #                    print 'data = ',data
#                        print "data['phot_g_mean_mag'] = ",data['phot_g_mean_mag']
                        dataList.append(np.float32(data['phot_g_mean_mag']))
                print 'len(dataList) = ',len(dataList)
                print 'filename ',filename,' read'
#print 'dataList = ',dataList
#print 'type(dataList) = ',type(dataList)
#print 'type(dataList[0]) = ',type(dataList[0])
print 'np.min(dataList) = ',np.min(dataList)
print 'np.max(dataList) = ',np.max(dataList)
#dataList = []
#for line in fileLines:
#    line_split = line.split(',')
#    data = {}
#    for iCol in range(len(header)):
#        data[header[iCol]] = line[iCol]
#    print 'data = ',data
#    dataList.append(data)
#print 'len(dataList) = ',len(dataList)

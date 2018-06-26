from collections import deque
import numpy as np
import os
import resource
import timeit

resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1))

#dataDir = '/Volumes/external/azuri/data/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/csv/'
dataDir = '/Users/azuri/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/csv/'
fileNameRoot = 'GaiaSource_%03d-%03d-%03d.csv'#% (release, tract, patch)
outDir = '/Volumes/external/azuri/data/gaia/lon-lat_temp/'
outFileRoot = 'GaiaSource_%d-%d_%d-%d.csv'# % (int(minLongitude), int(maxLongitude), int(minLatitude), int(maxLatitude))

filename = os.path.join(dataDir, fileNameRoot % (0, 0, 0))
start_time = timeit.default_timer()
with open(filename, 'rb') as f:
    for line in f:
        header = line.split(',')
        break
print 'header read in ', timeit.default_timer() - start_time,' s'
print 'header = ',header
lCol = -1
bCol = -1
for iHead in range(len(header)):
    if header[iHead] == 'l':
        lCol = iHead
    if header[iHead] == 'b':
        bCol = iHead
print 'lCol = ',lCol,', bCol = ',bCol
if lCol == -1:
    print 'lCol not found'
    STOP
if bCol == -1:
    print 'bCol not found'
    STOP

stepLong = 10.0
stepLat = 10.0
longitudes = np.arange(0.0, 361.0, stepLong)
print 'longitudes = ',longitudes
latitudes = np.arange(-90.0, 91.0, stepLat)
print 'latitudes = ',latitudes

nLongitudes = len(longitudes)
nLatitudes = len(latitudes)

outFiles = []
for iLong in np.arange(1, nLongitudes):
    files = []
    for iLat in np.arange(1, nLatitudes):
        outFileName = os.path.join(outDir,
                                   outFileRoot % (int(longitudes[iLong-1]),
                                                  int(longitudes[iLong]),
                                                  int(latitudes[iLat-1]),
                                                  int(latitudes[iLat])))
        file = open(outFileName, 'w')
        files.append(file)
    outFiles.append(files)
print 'len(outFiles) = ',len(outFiles)
print 'len(outFiles[0]) = ',len(outFiles[0])

for iRelease in range(3):
    for iTract in np.arange(0, 1000):
        for iPatch in np.arange(0, 1000):
            #start timer
            start_time = timeit.default_timer()

            filename = os.path.join(dataDir, fileNameRoot % (iRelease, iTract, iPatch))
            if os.path.exists(filename):
                with open(filename, 'rb') as f:
                    print 'filename ',filename,' opened in ', timeit.default_timer() - start_time,' s'
                    iLine = 0
    #                    file_content = f.read()
    #                    file_lines = file_content.splitlines()
    #                    del file_lines[0]
    #                    for line in file_lines:
                    for line in f:
    #                        print 'line ',iLine,' read at ', timeit.default_timer() - start_time,' s'
                        if iLine == 0:
                            iLine = 1
                            continue
#                            pass
                        line_split = line.split(',')

                        lon = np.float32(line_split[lCol])
                        lat = np.float32(line_split[bCol])
                        found = False
                        for iLong in np.arange(1, nLongitudes):
                            for iLat in np.arange(1, nLatitudes):
                                if (lon >= int(longitudes[iLong-1])
                                    and lon < int(longitudes[iLong])):
                                    #print 'longitude ',line_split[lCol],' found'
                                    if (lat >= int(latitudes[iLat-1])
                                        and lat < int(latitudes[iLat])):
                                        #print 'latitude ',line_split[bCol],' found'
                                        outFiles[iLong-1][iLat-1].write("%s\n" % line)
#                                        print 'writing lon = ',lon,', lat = ',lat,' to file outFiles[',iLong-1,'][',iLat-1,'] = ',outFiles[iLong-1][iLat-1].name
                                        found = True
                        if not found:
                            print 'ERROR: longitude ',line_split[lCol],', latitude ',line_split[bCol],' not found'
                            STOP
#                        if iLine == 10000:
#                            print 'read 10,000 lines in ', timeit.default_timer() - start_time,' s'
#                        iLine+=1
                print 'filename ',filename,' read in ', timeit.default_timer() - start_time,' s'

for files in outFiles:
    for file in files:
        file.close()
print 'all done'

from datetime import date, timedelta, datetime

def multikeysort(items, columns):
    from operator import itemgetter
    comparers = [((itemgetter(col[1:].strip()), -1) if col.startswith('-') else
                  (itemgetter(col.strip()), 1)) for col in columns]
    def comparer(left, right):
        for fn, mult in comparers:
            result = cmp(fn(left), fn(right))
            if result:
                return mult * result
        else:
            return 0

# --- fit y = a*x + b
# --- return value w: w[0] = a, w[1] = b
def linReg(x,y):
    A = np.array([ x, np.ones(len(x))])
#    print 'linReg: x.shape = ',x.shape,', A.shape = ',A.shape,', y.shape = ',y.shape
    w = []
    for iy in range(y.shape[0]):
        w.append(np.linalg.lstsq(A.T,y[iy,:])[0])
#        print 'linReg: iy = ',iy,': w[',iy,'] = ',w[iy]
    return w

# axis: 0 (columns) or 1 (rows)
# width: odd number
def boxCarMedianSmooth(imageData, axis, width):
    newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0,0]))
    if axis == 0:
        for iRow in range(imageData.shape[0]):
            for iCol in range(imageData.shape[1]):
                if iCol < int(width/2.0):
                    iColStart = 0
                    iColEnd = iCol+int(width/2.0)+1
                elif iCol > imageData.shape[1]-int(width/2.0)-1:
                    iColStart = iCol - int(width/2.0)
                    iColEnd = imageData.shape[1]
                else:
                    iColStart = iCol - int(width/2.0)
                    iColEnd = iCol + int(width/2.0) + 1
                newDataArray[iRow,iCol] = np.median(imageData[iRow,iColStart:iColEnd])
#                print 'iRow = ',iRow,', iCol = ',iCol,': iColStart = ',iColStart,', iColEnd = ',iColEnd,': imageData[iRow,iColStart:iColEnd] = ',imageData[iRow,iColStart:iColEnd],': median = ',newDataArray[iRow,iCol]
    elif axis == 1:
        for iCol in range(imageData.shape[1]):
            for iRow in range(imageData.shape[0]):
                if iRow < int(width/2.0):
                    iRowStart = 0
                    iRowEnd = iRow+int(width/2.0)+1
                elif iRow > imageData.shape[0]-int(width/2.0)-1:
                    iRowStart = iRow - int(width/2.0)
                    iRowEnd = imageData.shape[1]
                else:
                    iRowStart = iRow - int(width/2.0)
                    iRowEnd = iRow + int(width/2.0) + 1
                newDataArray[iRow,iCol] = np.median(imageData[iRowStart:iRowEnd,iCol])
#                print 'iCol = ',iCol,', iRow = ',iRow,': iRowStart = ',iRowStart,', iRowEnd = ',iRowEnd,': imageData[iRowStart:iRowEnd,iCol] = ',imageData[iRowStart:iRowEnd,iCol],': median = ',newDataArray[iRow,iCol]
    else:
        print 'ERROR: axis(=',axis,') out of bounds [0,1]'
    return newDataArray

# --- get date from string of the form yyyy-mm-ddThh:mm:ss.mmm
def getDate(dateStr):
    dateStr=dateStr[:dateStr.find('T')]
    day=date(int(dateStr[:dateStr.find('-')]),
             int(dateStr[dateStr.find('-')+1:dateStr.rfind('-')]),
             int(dateStr[dateStr.rfind('-')+1:]))
    print 'day = ',day
    return day

def getDateTime(dateTimeStr):
    dateStr = dateTimeStr[:dateTimeStr.find('T')]
    timeStr = dateTimeStr[dateTimeStr.find('T')+1:]
    time=datetime(int(dateStr[:dateStr.find('-')]),
             int(dateStr[dateStr.find('-')+1:dateStr.rfind('-')]),
             int(dateStr[dateStr.rfind('-')+1:]),
             int(timeStr[0:timeStr.find(':')]),
             int(timeStr[timeStr.find(':')+1:timeStr.rfind(':')]),
             int(timeStr[timeStr.rfind(':')+1:timeStr.find('.')]))
    print 'time = ',time
    return time

def findClosestInTime(time, times):
    diffMin = 100000000.
    timeOut = times[0]
    timeIdx = 0
    idx = 0
    for timeA in times:
        diff = time - timeA
        timediff = abs(diff.total_seconds())
        print 'diff = ',diff,', timediff = ',type(timediff),': ',timediff
        if timediff < diffMin:
            diffMin = timediff
            timeOut = timeA
            timeIdx = idx
            idx += 1
    return [timeOut, timeIdx, diffMin]
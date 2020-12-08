import numpy as np
import xlrd
import matplotlib.pyplot as plt

fNameList = '/Users/azuri/daten/paola/allFiles.txt'
outName = '/Users/azuri/daten/paola/results_%s.txt'

def boxCarMedianSmooth(inputData, boxSize):
    print('inputData.shape = ',inputData.shape)
    outputData = np.zeros(inputData.shape[0])
    print('outputData.shape = ',outputData.shape)
    outputData[0:int(boxSize/2)] = inputData[0:int(boxSize/2)]
    print("outputData[0:int(boxSize/2)] = ",outputData[0:int(boxSize/2)])
    for i in np.arange(int(boxSize/2),inputData.shape[0]-int(boxSize/2),1):
        outputData[i] = np.median(inputData[i-int(boxSize/2):i+int(boxSize/2)+1])
        print("outputData[",i,"] = np.median(inputData[i-int(boxSize/2):i+int(boxSize/2)+1]=",np.sort(inputData[i-int(boxSize/2):i+int(boxSize/2)+1]),"] = ",np.median(inputData[i-int(boxSize/2):i+int(boxSize/2)+1]),")")
    outputData[inputData.shape[0]-int(boxSize/2):] = inputData[inputData.shape[0]-int(boxSize/2):]
    print("outputData[inputData.shape[0]-int(boxSize/2):]")
#    STOP
    return outputData

with open(fNameList,'r') as f:
    lines= f.readlines()
lines = [line.strip() for line in lines]
print('lines = ',lines)

for fName in lines:
    print('fName = <'+fName+'>')
    wb = xlrd.open_workbook(fName)
    print('dir(wb) = ',dir(wb))

    for sheetNo in np.arange(3,wb.nsheets,1):
        sheet = wb.sheet_by_index(sheetNo)
        print('dir(sheet) = ',dir(sheet))

        print(sheet.show_sheet_headers)

        print(' ')
        if sheet.ncols > 4 and sheet.nrows > 4:
            if sheet.cell_value(5, 5) != '':

                print(' ')
                for i in range(sheet.ncols):
                    print(sheet.cell_value(5, i))

                print(' ')
                newValues = []
                for i in np.arange(4,sheet.nrows,1):
                    val = sheet.cell_value(i, 5)
                    if val != '':
                        print('i = ',i,': val = ',val)
                        newValues.append(-float(val))
                print('newValues = ',newValues)

                plt.plot(newValues)
                newValues = np.array(newValues)
                newValuesSmoothed = boxCarMedianSmooth(newValues,15)
                for i in range(newValues.shape[0]):
                    print('newValues[',i,'] = ',newValues[i],', newValuesSmootned[',i,'] = ',newValuesSmoothed[i])
                #STOP
                plt.plot(newValuesSmoothed)

                maxima = []
                minima = []
                indices = [0]
                oldValue = newValuesSmoothed[0]
                currentValue = newValuesSmoothed[0]
                nextValue = None
                for i in np.arange(1,newValuesSmoothed.shape[0]-1,1):
                    if newValuesSmoothed[i] == currentValue:
                        indices.append(i)
                    else:
                        nextValue = newValuesSmoothed[i]
                        print('i = ',i,': oldValue = ',oldValue,', currentValue = ',currentValue,', nextValue = ',nextValue,': newValuesSmoothed[indices=',indices,'] = ',newValuesSmoothed[indices])
                        if (newValuesSmoothed[indices[0]] < oldValue) and (newValuesSmoothed[indices[0]] < nextValue):
                            minima.append([int(np.mean(indices)),newValuesSmoothed[indices[0]]])
                            print('oldValue = ',oldValue,', currentValue = ',currentValue,', newValue = ',nextValue,': minima = ',minima)
                #            STOP
                        elif (newValuesSmoothed[indices[0]] > oldValue) and (newValuesSmoothed[indices[0]] > nextValue):
                            maxima.append([int(np.mean(indices)),newValuesSmoothed[indices[0]]])
                            print('i = ',i,': newValuesSmoothed[indices[0]] = ',newValuesSmoothed[indices[0]],', newValuesSmoothed[i] = ',newValuesSmoothed[i],', oldValue = ',oldValue,', currentValue = ',currentValue,', nextValue = ',nextValue,': maxima = ',maxima)
                #            STOP
                #        else:
                        oldValue = newValuesSmoothed[indices[0]]
                        currentValue = newValuesSmoothed[i]
                        indices = [i]

                print('minima = ',minima)
                print('maxima = ',maxima)

                plt.scatter([m[0] for m in minima], [m[1] for m in minima])
                plt.scatter([m[0] for m in maxima], [m[1] for m in maxima])
                plt.show()

                with open(outName % (fName[fName.rfind('/')+1:fName.rfind('.')]+'_sheet'+str(sheetNo)),'w') as f:
                    for min in minima:
                        print('writing ',min)
                        f.write('%d %.2f\n' % (min[0],min[1]))
                    for max in maxima:
                        print('writing ',max)
                        f.write('%d %.2f\n' % (max[0],max[1]))

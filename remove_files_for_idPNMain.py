import os
import sys

if __name__ == '__main__':
    allFilesName = 'allFindingCharts.list'
    os.system('ls */*/ > '+allFilesName)
    with open(allFilesName,'r') as f:
        allFiles = f.readlines()

    allFullPaths = []
    path = ''
    for i in range(len(allFiles)):
        line = allFiles[i].strip()
        if line != '':
            if ':' in line:
                path = line[:line.rfind(':')]
            else:
                allFullPaths.append(os.path.join(path,line))
    for fName in allFullPaths:
        if ('alt' in fName) and ('hashID' in fName):
            idPNMain = fName[fName.rfind('_')+1:]
        elif ('moonDist' in fName):
            tmp = fName[:fName.rfind('_')]
            idPNMain = tmp[tmp.rfind('_')+1:]
        elif len(fName[fName.rfind('/')+1]) < 6:
            idPNMain = fName[fName.rfind('/')+1:]
        else:
            print('fName = ',fName)
            STOP
        if idPNMain == sys.argv[1]:
            print('removing file '+fName)
            if os.path.exists(fName):
                os.system('rm -r '+fName)
            if os.path.islink(fName):
#                print('still exists, trying again')
                os.unlink(fName)

import os
import sys

if __name__ == '__main__':
    allFilesName = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/allFindingCharts.list'
    os.system('ls /Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/*/*/*/ > '+allFilesName)
    os.system('ls /Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/*/*/ >> '+allFilesName)
    observedDir = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/observed'
    os.system('mkdir '+observedDir)
    with open(allFilesName,'r') as f:
        allFiles = f.readlines()

    allFullPaths = []
    path = ''
    for i in range(len(allFiles)):
        line = allFiles[i].strip()
        if line != '':
            if line[-1] == ':':
                path = line[:line.rfind(':')]
                #print('path = ',path)
#                STOP
            else:
                allFullPaths.append(os.path.join(path,line))
    for fName in allFullPaths:
        tmp = fName[fName.rfind('/')+1:]
        if tmp[tmp.rfind('.'):] == '.png':
            #print('tmp = ',tmp,' ends with .png')
#            if 'RAJ2000' in fName:
#                STOP
            if ('alt' in tmp) and ('hashID' in tmp):
                idPNMain = tmp[tmp.find('_hashID_')+8:]
                idPNMain = idPNMain[:idPNMain.find('_')]
                #print('fName = ',fName,': idPNMain = ',idPNMain)
                #STOP
            elif ('moonDist' in tmp) and (not 'J2000' in tmp):
                tmp = tmp[tmp.find('_')+1:]
                tmp = tmp[tmp.find('_')+1:]
                idPNMain = tmp[:tmp.find('_')]
                #print('fName = ',fName,': idPNMain = ',idPNMain)
                #STOP
            elif ('J2000' in tmp) or (tmp[:7] == 'RAJ2000'):
                #print('tmp = ',tmp,' has J2000 in it')
                #print('tmp = <'+tmp+'>')
                tmp = tmp[tmp.find('HASH-ID=')+8:]
                #print('tmp = <'+tmp+'>')
                idPNMain = tmp[:tmp.find('_')]
                #print('fName = ',fName,': idPNMain = ',idPNMain)
#                STOP
    #        elif len(tmp[tmp.rfind('/')+1]) < 6:
    #            idPNMain = tmp[tmp.rfind('/')+1:]
    #            print('fName = ',fName,': idPNMain = ',idPNMain)
    #            STOP
            else:
                idPNMain = '0'
    #            STOP
            if idPNMain == sys.argv[1]:
                print('removing file '+fName)
                if True:
                    if os.path.exists(fName):
                        print('removing '+fName)
                        os.system('mv '+fName+' '+observedDir)
                    if os.path.islink(fName):
        #                print('still exists, trying again')
                        print('removing '+fName)
                        os.unlink(fName)

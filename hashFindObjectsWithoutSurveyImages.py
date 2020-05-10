import numpy as np
import os

import csvFree,csvData

imageDirList = '/Users/azuri/daten/uni/HKU/HASH/%s_images.txt'
hashPNMain = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full.csv'

commandFileOut = '/Users/azuri/daten/uni/HKU/HASH/addSurveyImages'
if os.path.exists(commandFileOut):
  os.remove(commandFileOut)

surveyCSV = csvData.CSVData()
surveyCSV.header = ['Name', 'lMin', 'lMax', 'bMin', 'bMax']
surveyCSV.append(['IPHAS', '29.', '215.', '-5.', '5.'])
surveyCSV.append(['VVV', '350.', '360.', '-14.5', '9.5'])
surveyCSV.append(['VVV', '0.', '10.', '-14.5', '9.5'])
surveyCSV.append(['VVV', '230.', '350.', '-4.5', '4.5'])
surveyCSV.append(['VVV', '10.', '20.', '-4.5', '4.5'])

if __name__ == '__main__':
    for iSurvey in range(surveyCSV.size()):
        with open(imageDirList % surveyCSV.getData('Name',iSurvey).lower(),'r') as f:
            imDirLines = f.readlines()

        idsWithoutSurveyImages = []
        freshIDFound = False
        idPNMain = 0
        idsStr = ''
        for line in imDirLines:
            if line[0] == '/':
                lastChar = 0-len(surveyCSV.getData('Name',iSurvey))-3
                idPNMain = line[len('/data/fermenter/PNImages/'):lastChar]
#                print('idPNMain = ',idPNMain)
                freshIDFound = True
            else:
                if line != '\n':
                    freshIDFound = False
                else:
                    if freshIDFound:
                        idsWithoutSurveyImages.append(idPNMain)
#        print('survey name = ',surveyCSV.getData('Name',iSurvey),': idsWithoutSurveyImages = ',len(idsWithoutSurveyImages),': ',idsWithoutSurveyImages)

        hashPN = csvFree.readCSVFile(hashPNMain)
        nMissing = 0
        for i in range(hashPN.size()):
            hashID = hashPN.getData('idPNMain',i)

            if hashID in idsWithoutSurveyImages:
                l = float(hashPN.getData('Glon',i))
                b = float(hashPN.getData('Glat',i))
                if ((l >= float(surveyCSV.getData('lMin',iSurvey)))
                    and (l <= float(surveyCSV.getData('lMax',iSurvey)))
                    and (b >= float(surveyCSV.getData('bMin',iSurvey)))
                    and (b <= float(surveyCSV.getData('bMax',iSurvey)))):
                    if idsStr == '':
                        idsStr = hashID
                    else:
                        idsStr += ','+hashID
                    nMissing += 1

        print('survey name = ',surveyCSV.getData('Name',iSurvey),': idsStr = ',idsStr)
        print('nMissing = ',nMissing)

        with open(commandFileOut,'a') as fCom:
            fCom.write('hashpn fetch '+surveyCSV.getData('Name',iSurvey).lower()+' '+idsStr+' -w force\n')
            fCom.write('hashpn brew all '+idsStr+'\n')



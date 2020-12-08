import csv

hashPNMainName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_241120.csv'
hashtbCNamesName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_241120.csv'
hashtbUsrCommName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbUsrComm_241120.csv'
hashFitsFilesName = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_241120.csv'
hashFrenchAmateursListName = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FRAall_HASHoutput_doubles_removed.csv'

outListTemp = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FRAall_HASHoutput_doubles_removed_CNames_spectraSources_comments.csv'
outList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FRAall_HASHoutput_doubles_removed_nonFRA_spectrum_found.csv'
outListAllSpec = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FRAall_HASHoutput_doubles_removed_spectrum_found.csv'

def createTempList():
    with open(outListTemp,'w') as fOut:
        fOut.write('idPNMain,RAJ2000,DECJ2000,PNstat,common names,spectra references,comments\n')
        with open(hashFrenchAmateursListName,'r') as hashFrenchAmateursListFile:
            hashFrenchAmateursList = csv.DictReader(hashFrenchAmateursListFile)
            for hashFrenchAmateursListRow in hashFrenchAmateursList:
                idPNMain = hashFrenchAmateursListRow['pndb']
                fOut.write(idPNMain + ',')
                with open(hashPNMainName, 'r') as hashPNMainFile:
                    hashPNMain = csv.DictReader(hashPNMainFile)
                    for hashPNMainRow in hashPNMain:
                        if (hashPNMainRow['idPNMain'] == idPNMain):
                            print('found idPNMain '+idPNMain+' in hashPNMain')
                            fOut.write(hashPNMainRow['RAJ2000']+','+hashPNMainRow['DECJ2000']+','+hashPNMainRow['PNstat']+',')
                with open(hashtbCNamesName, 'r') as hashtbCNamesFile:
                    hashtbCNames = csv.DictReader(hashtbCNamesFile)
                    found = False
                    for hashtbCNamesRow in hashtbCNames:
                        if (hashtbCNamesRow['idPNMain'] == idPNMain):
                            print('found idPNMain '+idPNMain+' in hashtbCNames')
                            if not found:
                                fOut.write(hashtbCNamesRow['Name'])
                            else:
                                fOut.write(';'+hashtbCNamesRow['Name'])
                            found = True
                fOut.write(',')
                with open(hashFitsFilesName, 'r') as hashFitsFilesFile:
                    hashFitsFiles = csv.DictReader(hashFitsFilesFile)
                    found = False
                    for hashFitsFilesRow in hashFitsFiles:
                        if (hashFitsFilesRow['idPNMain'] == idPNMain):
                            print('found idPNMain '+idPNMain+' in hashFitsFiles')
                            if (hashFitsFilesRow['reference'] != 'FrenchAmateurs'):
                                print('no FrenchAmateurs spectrum found')
                            if not found:
                                fOut.write(hashFitsFilesRow['reference'])
                            else:
                                fOut.write(';'+hashFitsFilesRow['reference'])
                            found = True
                fOut.write(',')
                with open(hashtbUsrCommName, 'r') as hashtbUsrCommFile:
                    hashtbUsrComm = csv.DictReader(hashtbUsrCommFile)
                    found = False
                    for hashtbUsrCommRow in hashtbUsrComm:
                        if (hashtbUsrCommRow['idPNMain'] == idPNMain):
                            print('found idPNMain '+idPNMain+' in hashtbUsrComm')
                            if not found:
                                fOut.write('"'+hashtbUsrCommRow['comment'].replace('\n',';')+'"')
                            else:
                                fOut.write(';"'+hashtbUsrCommRow['comment'].replace('\n',';')+'"')
                            found = True
                fOut.write('\n')

def removeFRASpectra():
    with open(outList,'w') as fOut:
        with open(outListAllSpec,'w') as fOutAll:
            fOut.write('idPNMain,RAJ2000,DECJ2000,PNstat,common names,spectra references,comments\n')
            fOutAll.write('idPNMain,RAJ2000,DECJ2000,PNstat,common names,spectra references,comments\n')
            with open(outListTemp, 'r') as tempFile:
                tempList = csv.DictReader(tempFile)
                for tempListRow in tempList:
                    refs = tempListRow['spectra references']
                    if refs != "":
                        print('refs = ',type(refs),': ',refs)
                        fOutAll.write(tempListRow['idPNMain'] + ',' + tempListRow['RAJ2000'] + ',' + tempListRow[
                            'DECJ2000'] + ',' + tempListRow['PNstat'] + ',' + tempListRow['common names'] + ',' +
                                   tempListRow['spectra references'] + ',"' + tempListRow['comments'].replace('"','') + '"\n')
                        if "FrenchAmateurs" not in refs:
                            fOut.write(tempListRow['idPNMain']+','+tempListRow['RAJ2000']+','+tempListRow['DECJ2000']+
                                       ','+tempListRow['PNstat']+','+tempListRow['common names']+','+
                                       tempListRow['spectra references']+',"' + tempListRow['comments'].replace('"','') + '"\n')
                            print('idPNMain = ',tempListRow['idPNMain'],': "FrenchAmateurs" NOT found in refs')
                        else:
                            print('idPNMain = ',tempListRow['idPNMain'],': "FrenchAmateurs" found in refs')

if __name__ == "__main__":
    createTempList()
    removeFRASpectra()
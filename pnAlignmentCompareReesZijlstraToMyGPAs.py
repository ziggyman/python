import os
import csvFree,csvData

path = '/Users/azuri/daten/uni/HKU/PN alignment'
rzFile = os.path.join(path, 'Rees_Zijlstra_table_with_HASH-ID.csv')
hashFile = os.path.join(path, 'HASH_true_PNe+004.2-05.9+005.9-02.6.csv')
myFile = os.path.join(path, 'PN-alignments.csv')

rzData = csvFree.readCSVFile(rzFile)
hashData = csvFree.readCSVFile(hashFile)
myData = csvFree.readCSVFile(myFile)

for i in range(rzData.size()):
    hashID = rzData.getData('HASH ID',i)
    rzGPA = float(rzData.getData('GPA',i))
    found = False
    for j in range(hashData.size()):
        if hashID == hashData.getData('idPNMain',j):
            print('HASH ID = ',hashID,': rzGPA = ',rzGPA,', HASH GPA = ',hashData.getData('GPA',j))
            continue
    for j in range(myData.size()):
        if hashID == myData.getData('HASH ID',j):
            myGPA = float(myData.getData('GPA',j))
            if myGPA > 180.:
                myGPA -= 180.
            if myGPA < 0.:
                myGPA += 180.
            print('HASH ID = ',hashID,': rzGPA = ',rzGPA,', my GPA = ',myGPA,', difference = ',myGPA-rzGPA)

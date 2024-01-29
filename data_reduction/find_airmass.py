import csvFree,csvData
from drUtils import getHeaderValue,setHeaderValue

areasFileName = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/CentralStars/R_data/2008-05-14/areas.csv'

areas = csvFree.readCSVFile(areasFileName)

for i in range(areas.size()):
    fNameCS = areas.getData('fName',i)
    if fNameCS[0] != '#':
        airmass = getHeaderValue(fNameCS,'AIRMASS')
        print('airmass = ',airmass)
        if airmass == None:
            fName = fNameCS.replace('CentralStars/','')
            fName = fName.replace('-sky','Ec')
            setHeaderValue(fNameCS,'AIRMASS',getHeaderValue(fName,'AIRMASS'))
            print('new airmass = ',getHeaderValue(fNameCS,'AIRMASS'))

        ra = getHeaderValue(fNameCS,'RA')
        print('ra = ',ra)
        if ra == None:
            fName = fNameCS.replace('CentralStars/','')
            fName = fName.replace('-sky','Ec')
            setHeaderValue(fNameCS,'RA',getHeaderValue(fName,'RA'))
            print('new RA = ',getHeaderValue(fNameCS,'RA'))

        dec = getHeaderValue(fNameCS,'DEC')
        print('dec = ',dec)
        if dec == None:
            fName = fNameCS.replace('CentralStars/','')
            fName = fName.replace('-sky','Ec')
            setHeaderValue(fNameCS,'DEC',getHeaderValue(fName,'DEC'))
            print('new DEC = ',getHeaderValue(fNameCS,'DEC'))

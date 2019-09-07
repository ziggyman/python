from findTargets import readCSV

#fNameIn = '/Users/azuri/daten/uni/HKU/observing/targets_SSO_priority_good_new.csv'
fNameIn = '/Users/azuri/daten/uni/HKU/observing/targets_SSO_new_noDiamGiven_good.csv'
#fNameIn = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_good.csv'
fNameOut = fNameIn[0:fNameIn.rfind('.')]+'_RA_DEC.dat'

lines = readCSV(fNameIn)

outList = []
with open(fNameOut,'w') as f:
    for line in lines:
        f.write(line['RAJ2000']+' '+line['DECJ2000']+'\n')

from findTargets import readCSV,writeCSV

fileA = '/Users/azuri/daten/uni/HKU/observing/targets_SSO_priority_new_good.csv'
fileB = '/Users/azuri/daten/uni/HKU/observing/targets_SSO_priority_good.csv'

fileNameOut = '/Users/azuri/daten/uni/HKU/observing/targets_SSO_priority_good_new.csv'

targetsA = readCSV(fileA)
targetsB = readCSV(fileB)

goodTargets = []

for lineA in targetsA:
    for lineB in targetsB:
        if lineA['idPNMain'] == lineB['idPNMain']:
            goodTargets.append(lineA)

writeCSV(goodTargets, fileNameOut, 'DRAJ2000')

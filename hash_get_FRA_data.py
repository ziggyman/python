import csv

fraFileName = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/RefFrenchPN_HASH_839obj_01022021.csv'
hashSpecFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_280121.csv'
pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_280121.csv'

tableOutFRA = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/RefFrenchPN_HASH_839obj_fitsDataFRA.csv'
tableOutNonFRA = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/RefFrenchPN_HASH_839obj_fitsDataNonFRA.csv'

def getObservatoryName(obs):
    if obs == 'SA':
        observatory = 'Sutherland (South Africa)'
    elif obs == 'KO':
        observatory = 'Kermerrien Observatory (Porspoder France)'
    elif obs == 'CN':
        observatory = 'OCA Calern (France)'
    elif (obs == 'CO') or (obs == 'CR'):
        observatory = 'Cornillon (France)'
    elif obs == 'AQ':
        observatory = 'Observatoire du pic de Château Renard (AstroQueyras; France)'
    elif obs == 'WH':
        observatory = 'WHT La Palma (Spain)'
    elif obs == 'SM':
        observatory = 'OAN-SPM; BC (Mexico)'
    elif obs == 'GT':
        observatory = "Observatoire Régional des Métiers Provence - Alpes - Côte d'Azur (France)"
    elif obs == 'HP':
        observatory = 'Observatoire de Haute-Provence (France)'
    elif obs == 'KP':
        observatory = 'Kitt Peak National Observatory (USA)'
    elif obs == 'OS':
        observatory = 'Observatorio de Sierra Nevada (Spain)'
    elif obs == 'FL':
        observatory = 'Australian Astronomical Observatory (Australia)'
    elif obs == 'IN':
        observatory = 'Isaac Newton Telescope (INT) La Palma (Spain)'
    elif obs == 'MY':
        observatory = 'Myslín (Czechia)'
    elif obs == 'MO':
        observatory = 'Mirranook Armidale (Australia)'
    elif obs == 'PM':
        observatory = 'Pic du Midi (France)'
    else:
        print('observatory <'+obs+'> not know')
        STOP

csvOutFra = []
csvOutNonFra = []
fraData = csv.DictReader(open(fraFileName))
fraIDs = []
for fraLine in fraData:
    idPNMain = fraLine['idPNMain']
    fraIDs.append(idPNMain)

""" All spectra of French Amateurs PN candidates not observed by them """
if False:
    fraData = csv.DictReader(open(fraFileName))
    i=0
    for fraLine in fraData:
        idPNMain = fraLine['idPNMain']
        hashSpecs = csv.DictReader(open(hashSpecFileName))
        for hashSpec in hashSpecs:
    #        print('hashSpec.items = ',list(hashSpec.items()))
    #        print('hashSpec.keys = ',list(hashSpec.keys()))
            if hashSpec['idPNMain'] == idPNMain:
                obs = hashSpec['fileName']
                if obs != 'LDu1_sum.fits':
                    print('i=',i,': fileName = ',obs)
                    obs = obs[obs.rfind('_')+1:obs.rfind('_')+3]
                    observatory = getObservatoryName(obs)
                    pnMain = csv.DictReader(open(pnMainFile))
                    for mainRow in pnMain:
                        if mainRow['idPNMain'] == idPNMain:
                            dictSpec = {'idPNMain': idPNMain,
                                        'name': fraLine['French Lists Name'],
                                        'RAJ2000': mainRow['RAJ2000'],
                                        'DECJ2000': mainRow['DECJ2000'],
                                        'DRAJ2000': mainRow['DRAJ2000'],
                                        'DDECJ2000': mainRow['DDECJ2000'],
                                        'date': hashSpec['dateObs'],
                                        'status': mainRow['PNstat'],
                                        'telescope': hashSpec['telescope'],
                                        'observatory': observatory
                                        }
                            if hashSpec['reference'] != 'FrenchAmateurs':
                                csvOutNonFra.append(dictSpec)
        i += 1

""" all spectra taken by French Amateurs """
j = 0
hashSpecs = csv.DictReader(open(hashSpecFileName))
for hashSpec in hashSpecs:
    if hashSpec['reference'] == 'FrenchAmateurs':
        idPNMain = hashSpec['idPNMain']
        if idPNMain not in fraIDs:
            obs = hashSpec['fileName']
            print('j=',j,': fileName = ',obs)
            obs = obs[obs.rfind('_')+1:obs.rfind('_')+3]
            observatory = getObservatoryName(obs)
            pnMain = csv.DictReader(open(pnMainFile))
            for mainRow in pnMain:
                if mainRow['idPNMain'] == idPNMain:
                    dictSpec = {'idPNMain': idPNMain,
                                'name': hashSpec['object'],
                                'RAJ2000': mainRow['RAJ2000'],
                                'DECJ2000': mainRow['DECJ2000'],
                                'DRAJ2000': mainRow['DRAJ2000'],
                                'DDECJ2000': mainRow['DDECJ2000'],
                                'date': hashSpec['dateObs'],
                                'status': mainRow['PNstat'],
                                'telescope': hashSpec['telescope'],
                                'observatory': observatory
                                }
                    if hashSpec['reference'] == 'FrenchAmateurs':
                        csvOutFra.append(dictSpec)
            j += 1


fieldnames = csvOutFra[0].keys()
print('fieldnames = ',fieldnames)
with open(tableOutFRA, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in csvOutFra:
        writer.writerow(row)
with open(tableOutNonFRA, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in csvOutNonFra:
        writer.writerow(row)

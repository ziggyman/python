import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.stats import circstats
import pycircstat

import csvFree,csvData
from pnOrientationUtils import calcMoments,calcMean,rose_plot,vectorDiagram,plotHammerProjection,linearOrderDiagram,plotEPA
from myUtils import hmsToDeg,dmsToDeg,findClosestObjectTo

fullSampleHASHFileName = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3_withGPA.csv'
fullSampleGPAsFileName = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/PN-alignments_finished.csv'
postCEPNeFileName = '/Users/azuri/daten/uni/HKU/publications/PNorientations/post-CEbinaryPNe.list'

fullSampleHASH = csvFree.readCSVFile(fullSampleHASHFileName)
fullSampleGPAs = csvFree.readCSVFile(fullSampleGPAsFileName)
postCEPNe = csvFree.readCSVFile(postCEPNeFileName,'\t',False)

show = True

nFound = 0
nNotFound = 0
idxs = []
for i in range(postCEPNe.size()):
    ra = hmsToDeg(postCEPNe.getData('RA',i))
    dec = dmsToDeg(postCEPNe.getData('DEC',i))
    idx,dist = findClosestObjectTo(ra,dec,fullSampleHASH,'DRAJ2000','DDECJ2000')
    if dist < 100.:
        print('found [',ra,', ',dec,'] at position ',idx,' with distance ',dist)
        idxs.append(idx[0])
        nFound += 1
    else:
        print('did not find ',ra,', ',dec,' in full HASH sample with GPAs, name = ',postCEPNe.getData('PN',i))
        nNotFound += 1
STOP
print('found ',nFound,' postCSPNe in full sample and did not find ',nNotFound)
idPNMainsHASH = fullSampleHASH.getData('idPNMain',idxs)
print('idPNMains = ',idPNMainsHASH)
idxGPAs = []
for i in range(len(idPNMainsHASH)):
    idx = fullSampleGPAs.find('HASH ID',idPNMainsHASH[i])
    if idx[0] >= 0:
        for id in idx:
            idxGPAs.append(id)
    else:
        print('ERROR: did not find idPNMain ',idPNMainsHASH[i])
        STOP
for flags in [[1,2,3],[1,2],[1]]:
    GPAs = []
    for i in range(len(idxGPAs)):
        if int(fullSampleGPAs.getData('flag',idxGPAs[i])) in flags:
            GPAs.append({'gpa': float(fullSampleGPAs.getData('GPA',idxGPAs[i])),
                        'flag': int(fullSampleGPAs.getData('flag',idxGPAs[i])),
                        'l': float(fullSampleGPAs.getData('l',idxGPAs[i])),
                        'b': float(fullSampleGPAs.getData('b',idxGPAs[i])),
                        })

    for i in range(len(GPAs)):
        if GPAs[i]['gpa'] < 0.:
            GPAs[i]['gpa'] += 180.
        if GPAs[i]['gpa'] > 180.:
            GPAs[i]['gpa'] -= 180.

    for sample in ['full','withoutBulge']:
        print('sample = '+sample)
        if sample == 'withoutBulge':
            for i in np.arange(len(GPAs)-1,-1,-1):
                if (np.abs(GPAs[i]['l']) < 10.) and (np.abs(GPAs[i]['b']) < 10.):
                    print('removing GPAs[',i,'] = ',GPAs[i])
                    del(GPAs[i])
        print('GPAs = ',len(GPAs),': ',GPAs)
        gpaRad = np.array([np.radians(2.*g['gpa']) for g in GPAs])
        pRad = circstats.rayleightest(gpaRad)
        print('pRad = ',pRad)

        moments = calcMoments([g['gpa'] for g in GPAs])
        if moments[0][0] < 0.:
            moments[0][0] += 180.
        print('moments = ',moments)

        rosePlotMean = moments[0][0]
        rosePlotSigma = moments[1][0]

        fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        rose_plot(ax, [g['gpa'] for g in GPAs], bins=16, density=True, offset=0, lab_unit="degrees",
                    start_zero=True, fill=True, color='blue', max_count=None,
                    max_size=None, smooth=False, mean=rosePlotMean, sigma=rosePlotSigma)
        fig.tight_layout()
        suffix = '_flags'
        for flag in flags:
            suffix += '%d' % (flag)
        print('suffix = <'+suffix+'>')
        suffix += '_'+sample
        print('suffix = <'+suffix+'>')
        fig.savefig('/Users/azuri/daten/uni/HKU/publications/PNorientations/post-CEbinaryPNe_rosePlot_'+suffix+'.png', bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close(fig)
        print('rosePlotMean = ',rosePlotMean,', rosePlotSigma = ',rosePlotSigma)

        vectorDiagram(np.array([g['gpa'] for g in GPAs])*2., fNameOut='/Users/azuri/daten/uni/HKU/publications/PNorientations/post-CEbinaryPNe_vectorDiagram_'+suffix+'.png')
        linearOrderDiagram([g['gpa'] for g in GPAs], fNameOut='/Users/azuri/daten/uni/HKU/publications/PNorientations/post-CEbinaryPNe_linearOrderDiagram_'+suffix+'.png')

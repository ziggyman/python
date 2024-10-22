from astropy.coordinates import EarthLocation
import numpy as np
import os
from drUtils import getImageData,writeFits1D,getHeader,dispCor,getListOfFiles,readFileToArr,calcResponse,applySensFuncs
from scipy import interpolate
from matplotlib import pyplot as plt
import astropy.io.fits as pyfits

path = '/Users/azuri/spectra/SAAO_June2024/080624/'
obsList = [os.path.join(path,'SCIENCE_HaTr5_a5671077_otzxfif.fits'),
           os.path.join(path,'SCIENCE_HaTr5_a5671078_otzxfif_cleaned.fits'),
           os.path.join(path,'SCIENCE_HaTr5_a5671079_otzxfif.fits'),
           ]
skyName = os.path.join(path,'SCIENCE_HaTr5-Sky_a5671215_otzxfif_cleaned.fits')
sky = getImageData(skyName,0)
print('sky.shape = ',sky.shape)
to_remove = [np.arange(798,809,1),np.arange(846,868,1),np.arange(927,931,1),np.arange(949,954,1)]

imSum = np.zeros(sky.shape)
for obsName in obsList:
    image = getImageData(obsName,0)
    imageOrig = getImageData(obsName,0)
#    z = np.ma.array(image, mask=mask)
#    x, y = np.mgrid[0:z.shape[0], 0:z.shape[1]]
#    x1 = x[~z.mask]
#    y1 = y[~z.mask]
#    z1 = z[~z.mask]
#    image = interpolate.interp2d(x1, y1, z1)(np.arange(z.shape[0]), np.arange(z.shape[1]))
#    xRange = np.arange(0,image.shape[1],1)
#    yRange = np.arange(0,image.shape[0],1)
#    idx = np.where((xRange >= np.min([cleanRange[0][0],cleanRange[1][0]])) & (xRange <= np.max([cleanRange[0][0],cleanRange[1][0]])))[0]
##                    print('clean: idx = ',idx)
#    idy = np.where((yRange >= np.min([cleanRange[0][1],cleanRange[1][1]])) & (yRange <= np.max([cleanRange[0][1],cleanRange[1][1]])))[0]
##                    print('clean: idy = ',idy)
#    for x in idx:
#        for y in idy:
#            image[y,x] = np.nan
##                            print('image[',y,',',x,'] = ',image[y,x])
##                    image[idx,idy] = np.nan
#    ax2DRect.cla()
##                    ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
##                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
    for i in range(image.shape[1]):
        vec = image[:,i]
        x = np.arange(0,image.shape[0],1)
        data = np.column_stack((x,vec))
        mask = np.ones(vec.shape[0], dtype=bool)
        remove = np.arange(26,31,1)
        mask[remove] = False
        newVec = data[mask,:]
        f = interpolate.interp1d(newVec[:,0],newVec[:,1])
        predicted = f(remove)
        vec[~mask] = predicted
        image[:,i] = vec
        imageOrig[:,i] = vec
    for remove in to_remove:
        for i in range(image.shape[0]):
            vec = image[i,:]
            x = np.arange(0,image.shape[1],1)
            data = np.column_stack((x,vec))
            mask = np.ones(vec.shape[0], dtype=bool)
            mask[remove] = False
            newVec = data[mask,:]
            f = interpolate.interp1d(newVec[:,0],newVec[:,1])
            predicted = f(remove)
            vec[~mask] = predicted
            image[i,:] = vec
    plt.imshow(image)
    plt.show()

    newSky = np.zeros(sky.shape)
    newSky[:,0:2] = sky[:,0:2] * np.mean(image[:,0:3]) / np.mean(sky[:,0:3])

    for i in np.arange(2,795):
        newSky[:,i] = sky[:,i] * np.mean(image[:,np.arange(i-2,i+3,1)]) / np.mean(sky[:,np.arange(i-2,i+3,1)])

    goodCols = []
    for i in np.arange(810,846,1):
        goodCols.append(i)
    for i in np.arange(870,924,1):
        goodCols.append(i)
    for i in np.arange(932,947,1):
        goodCols.append(i)
    newSky[:,795:957] = sky[:,795:957] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    for i in np.arange(957,1397,1):
        newSky[:,i] = sky[:,i] * np.mean(image[:,i-2:i+3]) / np.mean(sky[:,i-2:i+3])

    for i in np.arange(1459,sky.shape[1]-2,1):
        newSky[:,i] = sky[:,i] * np.mean(image[:,i-2:i+3]) / np.mean(sky[:,i-2:i+3])

    if False:
        goodCols = []
        for i in np.arange(1073,1093,1):
            goodCols.append(i)
        for i in np.arange(1107,1127,1):
            goodCols.append(i)
        #for i in np.arange(1073,1127,1):
        newSky[:,1073:1127] = sky[:,1073:1127] * np.mean(image[:,1073:1127]) / np.mean(sky[:,1073:1127]) #np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    goodCols = []
    for i in np.arange(1116,1136,1):
        goodCols.append(i)
    for i in np.arange(1150,1170,1):
        goodCols.append(i)
    #for i in np.arange(1116,1170,1):
    newSky[:,1116:1170] = sky[:,1116:1170] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])


    if False:
        goodCols = []
        for i in np.arange(1178,1198,1):
            goodCols.append(i)
        for i in np.arange(1212,1232,1):
            goodCols.append(i)
        #for i in np.arange(1178,1232,1):
        newSky[:,1178:1232] = sky[:,1178:1232] * np.median(image[:,1178:1232]) / np.median(sky[:,1178:1232])#np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])


    goodCols = []
    for i in np.arange(1308,1328,1):
        goodCols.append(i)
    for i in np.arange(1342,1362,1):
        goodCols.append(i)
    #for i in np.arange(1308,1362,1):
    newSky[:,1308:1362] = sky[:,1308:1362] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    goodCols = []
    for i in np.arange(1377,1395,1):
        goodCols.append(i)
    for i in np.arange(1406,1426,1):
        goodCols.append(i)
    #for i in np.arange(1377,1426,1):
    newSky[:,1377:1426] = sky[:,1377:1426] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

#    goodCols = []
#    for i in np.arange(1407,1445,1):
#        goodCols.append(i)
##    #for i in np.arange(1397,1459,1):
#    newSky[:,1397:1459] = sky[:,1397:1459] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    goodCols = []
    for i in np.arange(1425,1445,1):
        goodCols.append(i)
    for i in np.arange(1460,1480,1):
        goodCols.append(i)
    #for i in np.arange(1425,1480,1):
    newSky[:,1425:1480] = sky[:,1425:1480] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    goodCols = []
    for i in np.arange(1539,1554,1):
        goodCols.append(i)
    for i in np.arange(1569,1585,1):
        goodCols.append(i)
    #for i in np.arange(1539,1585,1):
    newSky[:,1539:1585] = sky[:,1539:1585] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])


    goodCols = []
    for i in np.arange(1604,1622,1):
        goodCols.append(i)
    for i in np.arange(1636,1653,1):
        goodCols.append(i)
    #for i in np.arange(1604,1653,1):
    newSky[:,1604:1653] = sky[:,1604:1653] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])

    goodCols = []
    for i in np.arange(1690,1704,1):
        goodCols.append(i)
    for i in np.arange(1718,1733,1):
        goodCols.append(i)
    #for i in np.arange(1690,1733,1):
    newSky[:,1690:1733] = sky[:,1690:1733] * np.mean(image[:,goodCols]) / np.mean(sky[:,goodCols])


    newSky[:,sky.shape[1]-2:sky.shape[1]] = sky[:,sky.shape[1]-2:sky.shape[1]] * np.mean(image[:,sky.shape[1]-3:sky.shape[1]]) / np.mean(sky[:,sky.shape[1]-3:sky.shape[1]])



#    skyRanges = [np.arange(0,76,1),
#                np.arange(76,137,1),
#                np.arange(137,1082,1),
#                np.arange(1082,1104,1),
#                np.arange(1104,1191,1),
#                np.arange(1191,1213,1),
#                np.arange(1213,image.shape[1],1)]
#    for skyRange in skyRanges:
#        sky[:,skyRange] = sky[:,skyRange] * np.mean(image[:,skyRange]) / np.mean(sky[:,skyRange])
    imageOrig -= newSky
    plt.imshow(imageOrig)
    plt.show()

    hdulist = pyfits.open(obsName)
    hdulist[0].data = imageOrig
    imMinusSkyName = obsName[:obsName.rfind('.')]+"-newSkyIm.fits"
    hdulist.writeto(imMinusSkyName,overwrite=True)
    hdulist.close()

    spec = np.zeros(imageOrig.shape[1])
    for i in range(spec.shape[0]):
        spec[i] = np.sum(imageOrig[:,i])
    writeFits1D(spec,imMinusSkyName[:imMinusSkyName.rfind('.')]+'Ec.fits',wavelength=np.arange(1,len(spec)+1,1),header=getHeader(obsName),CRVAL1=None,CRPIX1=None,CDELT1=None)
    plt.plot(spec)
    plt.show()
    imSum += imageOrig

imSum = imSum / 3.
hdulist = pyfits.open(obsList[1])
hdulist[0].data = imageOrig
imMinusSkyName = obsName[:obsName.rfind('/')]+"/HaTr5_SA080624_combined.fits"
hdulist.writeto(imMinusSkyName,overwrite=True)
hdulist.close()

spec = np.zeros(imSum.shape[1])
for i in range(spec.shape[0]):
    spec[i] = np.sum(imSum[:,i])
writeFits1D(spec,imMinusSkyName[:imMinusSkyName.rfind('/')]+'/HaTr5_SA080624_combinedEc.fits',wavelength=np.arange(1,len(spec)+1,1),header=getHeader(obsName),CRVAL1=None,CRPIX1=None,CDELT1=None)
plt.plot(spec)
plt.show()

inputList = getListOfFiles(os.path.join(path,'arc_otzxfifEc.list'))
wavelengthsOrig = []
for i in range(len(inputList)):
    wLenStr = readFileToArr(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat')
    wLens = [float(wLen) for wLen in wLenStr]
    wavelengthsOrig.append(np.asarray(wLens))

plt.plot(wavelengthsOrig[14],spec)
plt.title('with wavelengthsOrig[14]')
plt.show()

observatoryLocation = EarthLocation.of_site('SAAO')
dispCor([os.path.join(path,imMinusSkyName[:imMinusSkyName.rfind('/')]+'/HaTr5_SA080624_combinedEc.fits')],
        getListOfFiles(os.path.join(path,'arc_otzxfiEc.list')),
        wavelengthsOrig,
        [os.path.join(path,imMinusSkyName[:imMinusSkyName.rfind('/')]+'/HaTr5_SA080624_combinedEcd.fits')],
        observatoryLocation,
        'TELRA',#'RA',#
        'TELDEC',#'DEC',#
        'DATE-OBS',
        doHelioCor = True)

sensFuncs, stdsDLam = calcResponse(os.path.join(path,'fluxstds_otzxfif.list'),
                        #getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
                        #wavelengthsOrig,
                        #areas,
                        stdStarNameEndsBefore='_a',
                        display=True)
print('sensFuncs = ',sensFuncs)
with open(os.path.join(path,'fluxstds_otzxfif.list'),'r') as f:
    fluxstds = f.readlines()
print('fluxstds = ',fluxstds)
fluxstds = [fl.strip() for fl in fluxstds]
print('fluxstds = ',fluxstds)
with open(os.path.join(path,'stdsDLambda.dat'),'w') as fl:
    for iStd in range(len(sensFuncs)):
        print('type(sensFuncs[iStd]) = ',type(sensFuncs[iStd]))
        print('sensFuncs[',iStd,'] = ',sensFuncs[iStd])
        sensFuncs[iStd].write(fluxstds[iStd][:fluxstds[iStd].rfind('.')]+'_sens.ecsv',overwrite=True)
        #sens = ascii.read(fluxstds[iStd][:fluxstds[iStd].rfind('.')]+'_sens.ecsv')
        #print('sens = ',sens)
        fl.write('%.3f\n' % (stdsDLam[iStd]))
        plt.plot(sensFuncs[iStd]['wave'],sensFuncs[iStd]['S'],label=fluxstds[iStd][fluxstds[iStd].rfind('/')+1:fluxstds[iStd].rfind('.')])
    plt.legend()
    plt.show()
if True:
    applySensFuncs([os.path.join(path,imMinusSkyName[:imMinusSkyName.rfind('/')]+'/HaTr5_SA080624_combinedEcd.fits')],
                    [os.path.join(path,imMinusSkyName[:imMinusSkyName.rfind('/')]+'/HaTr5_SA080624_combinedEcdF.fits')],
                    stdsDLam,
                    sensFuncs,
                    readFileToArr(os.path.join(path,'fluxstds_otzxfif.list')))

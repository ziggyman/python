import numpy as np
from matplotlib import pyplot as plt
from drUtils import getImageData,getWavelengthArr,getHeader,getHeaderValue,load_wcs_from_file
from myUtils import getArcsecDistance,getXYFromRaDec,getRaDecFromXY,angularDistancePyAsl,hmsToDeg,dmsToDeg,applyVRadCorrection

saveFigs = False
vrad = -120.

specFileBlue = '/Users/azuri/daten/uni/HKU/Kamila/OIII_ratios/Red and blue cube/Blue_final_cube_angstroms.fits'
specFileRed = '/Users/azuri/daten/uni/HKU/Kamila/OIII_ratios/Red and blue cube/Red_final_cube_angstroms.fits'

#blue
header = getHeader(specFileBlue,1)
print('header[CTYPE1] = ',header['CTYPE1'])
print('header[CTYPE2] = ',header['CTYPE2'])
print('header[CTYPE3] = ',header['CTYPE3'])

centerBlue = getXYFromRaDec(specFileBlue,'08:46:53.50','17:52:45.50')#[16.29,25.7]
print('centerBlue = ',centerBlue)
xCenterBlue = 14.9
yCenterBlue = 26.01
xDistBlue = 8.91
yDistBlue = 27.01
dist = getArcsecDistance(specFileBlue,xCenterBlue,yCenterBlue,xDistBlue,yDistBlue)
print('dist = ',dist)

raBlue,decBlue = getRaDecFromXY(specFileBlue,xCenterBlue,yCenterBlue)
print('raBlue = ',raBlue,', decBlue = ',decBlue)
raBlueDeg = hmsToDeg(raBlue)
decBlueDeg = dmsToDeg(decBlue)
print('raBlueDeg = ',raBlueDeg)
print('decBlueDeg = ',decBlueDeg)

nxBlue = getHeaderValue(specFileBlue,'NAXIS1',1)
nyBlue = getHeaderValue(specFileBlue,'NAXIS2',1)
nwBlue = getHeaderValue(specFileBlue,'NAXIS3',1)

blueSpec2D = getImageData(specFileBlue,1)
print('blueSpec2D.shape = ',blueSpec2D.shape)

#red
xCenterRed = 16.17
yCenterRed = 25.16

#xCenterRed,yCenterRed = getXYFromRaDec(specFileRed,raRed,decRed)#[14.5,25.88]
print('centerRed = [',xCenterRed,',',yCenterRed,']')
#raRed,decRed = getRaDecFromXY(specFileRed,14.5,25.88)
#print('raRed = ',raRed,', decRed = ',decRed)
redSpec2D = getImageData(specFileRed,1)
dist = 0.


lines = [[3874.,'[Ne III]'],
         [3969.,'[Ne III]'],
         [4026.,'He I'],
         [4072.,'[S II]'],
         [4102.,'H_delta'],
         [4269.,'C II'],
         [4340.,'H_gamma'],
         [4363.,'[O III]'],
         [4388.,'He I'],
         [4472.,'He I'],
         [4542.,'He II'],
         [4571.,'[Mg I]'],
         [4648.,'C III'],
         [4686.,'He II'],
         [4723.,'[Ne IV]'],
         [4740.,'[Ar IV]'],
         [4861.,'H_beta'],
         [4922.,'He I'],
         [4959.,'[O III]'],
         [5007.,'[O III]'],
         [5199.,'N I'],
         [5412.,'He II'],
         [5518.,'[Cl III]'],
         [5538.,'[Cl III]'],
         [5577.,'[O I]'],
         [5754.,'[N II]'],
         [5876.,'He I'],
         [6300.,'[O I]'],
         [6312.,'[S III]'],
         [6364.,'[O I]'],
         [6435.,'[Ar V]'],
         [6548.,'[N II]'],
         [6563.,'H_alpha'],
         [6583.,'[N II]'],
        ]

for arm in ['blue','red']:
    if arm == 'blue':
        specFile = specFileBlue
        armSpec2D = blueSpec2D
        xCenter = xCenterBlue
        yCenter = yCenterBlue
        wcs = load_wcs_from_file(specFile,1)
        world = wcs.wcs_pix2world(xCenter,yCenter,0,0)
        print('world = ',world)
        worldRA = world[0]
        worldDEC = world[1]
        print('worldRA = ',worldRA)
        print('worldDEC = ',worldDEC)

        worldDist = wcs.wcs_pix2world(xDistBlue,yDistBlue,0,0)
        print('worldDist = ',worldDist)
        worldDistRA = worldDist[0]
        worldDistDEC = worldDist[1]
        print('worldDistRA = ',worldDistRA)
        print('worldDistDEC = ',worldDistDEC)

        dist = angularDistancePyAsl(worldRA,worldDEC,worldDistRA,worldDistDEC)*3600.
        print('dist = ',dist)
    else:
        print('dist = ',dist)
        #STOP
        specFile = specFileRed
        armSpec2D = redSpec2D
        xCenter = xCenterRed
        yCenter = yCenterRed


    wLen = getWavelengthArr(specFile,1,3)
    print('wLen = ',len(wLen),': ',wLen[0],' - ',wLen[len(wLen)-1])

    wcs = load_wcs_from_file(specFile,1)
    print('wcs = ',wcs)
    world = wcs.wcs_pix2world(xCenter,yCenter,0,0)
    print('world = ',world)
    worldRA = world[0]
    worldDEC = world[1]


    plt.imshow(armSpec2D[0,:,:])
    plt.plot([xCenter],[yCenter],'r+')
    if saveFigs:
        plt.savefig(specFile[:specFile.rfind('.')]+'_2D_wLen=%.2f.png' % (wLen[0]), bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    nx = getHeaderValue(specFile,'NAXIS1',1)
    ny = getHeaderValue(specFile,'NAXIS2',1)
    nw = getHeaderValue(specFile,'NAXIS3',1)

    print(arm,': nx = ',nx,', ny = ',ny,', nw = ',nw)

    world = wcs.wcs_pix2world(xCenter,yCenter,0,0)
    print('world = ',world)
    worldRA = world[0]
    worldDEC = world[1]
    print('worldRA = ',worldRA)
    print('worldDEC = ',worldDEC)

    wLen = getWavelengthArr(specFile,1,3)

    dist2D = np.zeros([ny,nx])

    for x in range(nx):
        for y in range(ny):
            ra,dec,lam = wcs.wcs_pix2world(x,y,0,0)
            dist2D[y,x] = angularDistancePyAsl(ra,dec,worldRA,worldDEC)*3600.
    plt.imshow(dist2D)
    if saveFigs:
        plt.savefig(specFile[:specFile.rfind('.')]+'_dists.png', bbox_inches='tight')
        plt.close()
    else:
        plt.show()

#    plt.imshow(armSpec2D[0,:,:])
#    plt.plot([xCenter],[yCenter],'b+')
#    for x in range(nx):
#        for y in range(ny):
#            if dist2D[y,x] <= dist:
#                plt.plot([x],[y],'r+')
#                print('[',x,',',y,'] is inside')
#    if saveFigs:
#        plt.savefig(specFile[:specFile.rfind('.')]+'_2D_wLen=%.2f_CS_marked.png' % (wLen[0]), bbox_inches='tight')
#        plt.close()
#    else:
#        plt.show()

    spec = np.zeros(len(wLen))
    plt.imshow(armSpec2D[0,:,:])
    for i in range(len(wLen)):
        for x in range(nx):
            for y in range(ny):
                if dist2D[y,x] < dist:
                    spec[i] += armSpec2D[i,x,y]
                    if i == 0:
                        print('pixel [',x,',',y,'] is inside')
                        plt.plot([x],[y],'b+')
    if saveFigs:
        plt.savefig(specFile[:specFile.rfind('.')]+'_2D_wLen=%.2f_CS_marked.png' % (wLen[0]), bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    plt.plot(wLen,spec)
    plt.title('original')
    if saveFigs:
        plt.savefig(specFile[:specFile.rfind('.')]+'_spec_orig.png', bbox_inches='tight')
        plt.close()
    else:
        plt.show()


    wLen = applyVRadCorrection(wLen,vrad)
    plt.plot(wLen,spec)
    plt.title('corrected for vrad=120 km/s')
    if saveFigs:
        plt.savefig(specFile[:specFile.rfind('.')]+'_spec_corrected_for_vrad=%d.png' % (vrad), bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    #STOP

    for line in lines:
        if (line[0] > wLen[0]) and (line[0] < wLen[len(wLen)-1]):
            spec2D = np.zeros(armSpec2D[0,:,:].shape)
            for i in range(len(wLen)):
                if np.abs(wLen[i] - line[0]) < 3.:
                    spec2D += armSpec2D[i,:,:]
            spec2D[np.isnan(spec2D)] = 0.
            print('max = ',np.max(spec2D))
            plt.imshow(spec2D,vmin=0.,vmax=np.mean(spec2D)*2.,origin='lower')
            plt.title('%s %d' % (line[1],line[0]))
            if saveFigs:
                plt.savefig(specFile[:specFile.rfind('.')]+'_2D_%s_%d.png' % (line[1],line[0]), bbox_inches='tight')
                plt.close()
            else:
                plt.show()

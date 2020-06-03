exec(open("/Users/azuri/entwicklung/python/pyraf/pa30_plot_vrad.py").read())# import getDate, findClosestDate,...
import math
import astropy.io.fits as pyfits
from myUtils import getWavelength

speedOfLight = 299792.458# km/s
rotatorAngleGTC = -60.0#-81.8436672903833+185.0

tabGTC = tab
SII6716aGTC = SII6716a
SII6716bGTC = SII6716b
SII6731aGTC = SII6731a
SII6731bGTC = SII6731b
ArIII7136GTC = ArIII7136
lambda0_OIII5007 = 5006.843
lambda0_Halpha6563 = 6562.79
distGTC = dist

exec(open("/Users/azuri/entwicklung/python/pyraf/pa30_sparsepak_data.py").read())

#plt.plot([0.0, 0.0], [0.0, 9000], 'g')
distX = [0.0-a['centerDistanceX'] for a in tab]
distY = [a['centerDistanceY'] for a in tab]
SII6716a = [a['SII6716a'][0] for a in tab]
SII6731a = [a['SII6731a'][0] for a in tab]
SII6716b = [a['SII6716b'][0] for a in tab]
SII6731b = [a['SII6731b'][0] for a in tab]
OIII5007a = [a['OIII5007a'][0] for a in tab]
OIII5007b = [a['OIII5007b'][0] for a in tab]
Halpha6563a = [a['Halpha6563a'][0] for a in tab]
Halpha6563b = [a['Halpha6563b'][0] for a in tab]

SII6716a_plot = []
SII6716b_plot = []
SII6731a_plot = []
SII6731b_plot = []
OIII5007a_plot = []
OIII5007b_plot = []
Halpha6563a_plot = []
Halpha6563b_plot = []
distX_SII6716a_plot = []
distY_SII6716a_plot = []
distX_SII6716b_plot = []
distY_SII6716b_plot = []
distX_SII6731a_plot = []
distY_SII6731a_plot = []
distX_SII6731b_plot = []
distY_SII6731b_plot = []
distX_OIII5007a_plot = []
distY_OIII5007a_plot = []
distX_OIII5007b_plot = []
distY_OIII5007b_plot = []
distX_Halpha6563a_plot = []
distY_Halpha6563a_plot = []
distX_Halpha6563b_plot = []
distY_Halpha6563b_plot = []
for i in range(len(SII6716a)):
    if SII6716a[i] != 0.0:
        SII6716a_plot.append(SII6716a[i])
        distX_SII6716a_plot.append(distX[i])
        distY_SII6716a_plot.append(distY[i])
    if SII6716b[i] != 0.0:
        SII6716b_plot.append(SII6716b[i])
        distX_SII6716b_plot.append(distX[i])
        distY_SII6716b_plot.append(distY[i])
    if SII6731a[i] != 0.0:
        SII6731a_plot.append(SII6731a[i])
        distX_SII6731a_plot.append(distX[i])
        distY_SII6731a_plot.append(distY[i])
    if SII6731b[i] != 0.0:
        SII6731b_plot.append(SII6731b[i])
        distX_SII6731b_plot.append(distX[i])
        distY_SII6731b_plot.append(distY[i])
    if OIII5007a[i] != 0.0:
        print('i=',i,': OIII5007a[i] = ',OIII5007a[i])
        OIII5007a_plot.append(OIII5007a[i])
        distX_OIII5007a_plot.append(distX[i])
        distY_OIII5007a_plot.append(distY[i])
    if OIII5007b[i] != 0.0:
        print('i=',i,': OIII5007b[i] = ',OIII5007b[i])
        OIII5007b_plot.append(OIII5007b[i])
        distX_OIII5007b_plot.append(distX[i])
        distY_OIII5007b_plot.append(distY[i])
    if Halpha6563a[i] != 0.0:
        print('i=',i,': Halpha6563a[i] = ',Halpha6563a[i])
        Halpha6563a_plot.append(Halpha6563a[i])
        distX_Halpha6563a_plot.append(distX[i])
        distY_Halpha6563a_plot.append(distY[i])
    if Halpha6563b[i] != 0.0:
        print('i=',i,': Halpha6563b[i] = ',Halpha6563b[i])
        Halpha6563b_plot.append(Halpha6563b[i])
        distX_Halpha6563b_plot.append(distX[i])
        distY_Halpha6563b_plot.append(distY[i])
print('SII6731b_plot = ',SII6731b_plot)

vradSII6716a = (np.array(SII6716a_plot) - lambda0SII6716) * speedOfLight / lambda0SII6716
vradSII6716b = (np.array(SII6716b_plot) - lambda0SII6716) * speedOfLight / lambda0SII6716

vradSII6731a = (np.array(SII6731a_plot) - lambda0SII6731) * speedOfLight / lambda0SII6731
vradSII6731b = (np.array(SII6731b_plot) - lambda0SII6731) * speedOfLight / lambda0SII6731

vradOIII5007a = (np.array(OIII5007a_plot) - lambda0_OIII5007) * speedOfLight / lambda0_OIII5007
vradOIII5007b = (np.array(OIII5007b_plot) - lambda0_OIII5007) * speedOfLight / lambda0_OIII5007

vradHalpha6563a = (np.array(Halpha6563a_plot) - lambda0_Halpha6563) * speedOfLight / lambda0_Halpha6563
vradHalpha6563b = (np.array(Halpha6563b_plot) - lambda0_Halpha6563) * speedOfLight / lambda0_Halpha6563

print('vradOIII5007a = ',vradOIII5007a)
print('vradOIII5007b = ',vradOIII5007b)
print('vradHalpha6563a = ',vradHalpha6563a)
print('vradHalpha6563b = ',vradHalpha6563b)
print('min(vradSII6716a) = ',min(vradSII6716a))
print('min(vradSII6731a) = ',min(vradSII6731a))
print('min(vradOIII5007a) = ',min(vradOIII5007a))
print('min(vradSII6716b) = ',min(vradSII6716b))
print('min(vradSII6731b) = ',min(vradSII6731b))
print('min(vradOIII5007b) = ',min(vradOIII5007b))
minval = min([min(vradSII6731a), min(vradSII6716a), min(vradSII6716b), min(vradSII6731b), min(vradOIII5007a), min(vradOIII5007b), min(vradHalpha6563a), min(vradHalpha6563b)])
maxval = max([max(vradSII6731a), max(vradSII6716a), max(vradSII6716b), max(vradSII6731b), max(vradOIII5007a), max(vradOIII5007b), max(vradHalpha6563a), max(vradHalpha6563b)])
print('minval = ',minval,', maxval = ',maxval)

cm = plt.cm.get_cmap('rainbow')
sc = plt.scatter(distX_SII6716a_plot,
                 distY_SII6716a_plot,
                 c=vradSII6716a,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] 6716')
plt.colorbar(sc, label='radial velocity [km/s]')

print('plt.ylim = ',plt.ylim)
plt.ylim = [-80., 80.]
print('plt.ylim = ',plt.ylim)
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SII6716a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_SII6716b_plot,
                 distY_SII6716b_plot,
                 c=vradSII6716b,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] 6716')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SII6716b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_SII6731a_plot,
                 distY_SII6731a_plot,
                 c=vradSII6731a,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] 6731')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SII6731a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_SII6731b_plot,
                 distY_SII6731b_plot,
                 c=vradSII6731b,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] 6731')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SII6731b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



sc = plt.scatter(distX_OIII5007a_plot,
                 distY_OIII5007a_plot,
                 c=vradOIII5007a,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[OIII] 5007')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_OIII5007a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_OIII5007b_plot,
                 distY_OIII5007b_plot,
                 c=vradOIII5007b,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[OIII] 5007')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_OIII5007b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



sc = plt.scatter(distX_Halpha6563a_plot,
                 distY_Halpha6563a_plot,
                 c=vradHalpha6563a,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='H-alpha 6563')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_Halpha6563a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_Halpha6563b_plot,
                 distY_Halpha6563b_plot,
                 c=vradHalpha6563b,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='H-alpha 6563')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_Halpha6563b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()





SIIa_plot = []
SIIb_plot = []
distX_SIIa_plot = []
distY_SIIa_plot = []
distX_SIIb_plot = []
distY_SIIb_plot = []
vradSIIa = []
vradSIIb = []
for i in range(len(SII6716a)):
    if (SII6716a[i] != 0.0) and (SII6731a[i] != 0.0):
        vradSIIa.append((((SII6716a[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((SII6731a[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)) / 2.0)
        SIIa_plot.append((SII6716a[i] + SII6731a[i]) / 2.0)
        distX_SIIa_plot.append(distX[i])
        distY_SIIa_plot.append(distY[i])
    elif SII6716a[i] != 0.0:
        vradSIIa.append((SII6716a[i] - lambda0SII6716) * speedOfLight / lambda0SII6716)
        SIIa_plot.append(SII6716a[i])
        distX_SIIa_plot.append(distX[i])
        distY_SIIa_plot.append(distY[i])
    elif SII6731a[i] != 0.0:
        vradSIIa.append((SII6731a[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)
        SIIa_plot.append(SII6731a[i])
        distX_SIIa_plot.append(distX[i])
        distY_SIIa_plot.append(distY[i])
    if (SII6716b[i] != 0.0) and (SII6731b[i] != 0.0):
        vradSIIb.append((((SII6716b[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((SII6731b[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)) / 2.0)
        SIIb_plot.append((SII6716b[i] + SII6731b[i]) / 2.0)
        distX_SIIb_plot.append(distX[i])
        distY_SIIb_plot.append(distY[i])
    elif SII6716b[i] != 0.0:
        vradSIIb.append((SII6716b[i] - lambda0SII6716) * speedOfLight / lambda0SII6716)
        SIIb_plot.append(SII6716b[i])
        distX_SIIb_plot.append(distX[i])
        distY_SIIb_plot.append(distY[i])
    elif SII6731b[i] != 0.0:
        vradSIIb.append((SII6731b[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)
        SIIb_plot.append(SII6731b[i])
        distX_SIIb_plot.append(distX[i])
        distY_SIIb_plot.append(distY[i])

print('min(vradSIIa) = ',min(vradSIIa))
minval = min([min(vradSIIa), min(vradSIIb)])
maxval = max([max(vradSIIa), max(vradSIIb)])
print('minval = ',minval,', maxval = ',maxval)

cm = plt.cm.get_cmap('rainbow')
sc = plt.scatter(distX_SIIa_plot,
                 distY_SIIa_plot,
                 c=vradSIIa,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII]')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SIIa.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_SIIb_plot,
                 distY_SIIb_plot,
                 c=vradSIIb,
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII]')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_SIIb.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

vradaGTC = []
vradbGTC = []
distXaGTC = []
distYaGTC = []
distXbGTC = []
distYbGTC = []
for i in range(len(distGTC)):
    distXTemp = distGTC[i] * math.cos(math.radians(rotatorAngleGTC))
    distYTemp = distGTC[i] * math.sin(math.radians(rotatorAngleGTC))
    print('distGTC[',i,'] = ',distGTC[i],': distXTemp = ',distXTemp,', distYTemp = ',distYTemp)
    if (SII6716aGTC[i] != 0.0) and (SII6731aGTC[i] != 0.0) and (ArIII7136GTC[i] != 0.0):
        vradaGTC.append((((SII6716aGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((SII6731aGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731) + ((ArIII7136GTC[i] - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136)) / 3.0)
        print('all 3: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif (SII6716aGTC[i] != 0.0) and (SII6731aGTC[i] != 0.0):
        vradaGTC.append((((SII6716aGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((SII6731aGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)) / 2.0)
        print('SII6716 and SII6731: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif (SII6716aGTC[i] != 0.0) and (ArIII7136GTC[i] != 0.0):
        vradaGTC.append((((SII6716aGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((ArIII7136GTC[i] - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136)) / 2.0)
        print('SII6716 and ArII7136: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif (SII6731aGTC[i] != 0.0) and (ArIII7136GTC[i] != 0.0):
        vradaGTC.append((((SII6731aGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731) + ((ArIII7136GTC[i] - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136)) / 2.0)
        print('SII6731 and ArII7136: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif SII6716aGTC[i] != 0.0:
        vradaGTC.append((SII6716aGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716)
        print('SII6716: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif SII6731aGTC[i] != 0.0:
        vradaGTC.append((SII6731aGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)
        print('SII6731: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)
    elif ArIII7136GTC[i] != 0.0:
        vradaGTC.append((ArIII7136GTC[i] - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136)
        print('ArII7136: vradaGTC[',len(vradaGTC)-1,'] = ',vradaGTC[len(vradaGTC)-1])
        distXaGTC.append(distXTemp)
        distYaGTC.append(distYTemp)

    if (SII6716bGTC[i] != 0.0) and (SII6731bGTC[i] != 0.0):
        vradbGTC.append((((SII6716bGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716) + ((SII6731bGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)) / 2.0)
        print('SII6716 and SII6731: vradbGTC[',len(vradbGTC)-1,'] = ',vradbGTC[len(vradbGTC)-1])
        distXbGTC.append(distXTemp)
        distYbGTC.append(distYTemp)
    elif SII6716bGTC[i] != 0.0:
        vradbGTC.append((SII6716bGTC[i] - lambda0SII6716) * speedOfLight / lambda0SII6716)
        print('SII6716: vradbGTC[',len(vradbGTC)-1,'] = ',vradbGTC[len(vradbGTC)-1])
        distXbGTC.append(distXTemp)
        distYbGTC.append(distYTemp)
    elif SII6731bGTC[i] != 0.0:
        vradbGTC.append((SII6731bGTC[i] - lambda0SII6731) * speedOfLight / lambda0SII6731)
        print('SII6731: vradbGTC[',len(vradbGTC)-1,'] = ',vradbGTC[len(vradbGTC)-1])
        distXbGTC.append(distXTemp)
        distYbGTC.append(distYTemp)

print('vradaGTC = ',vradaGTC)
print('vradbGTC = ',vradbGTC)

minval = min([min(vradSIIa), min(vradSIIb), min(vradaGTC)])
maxval = max([max(vradSIIa), max(vradSIIb), max(vradbGTC)])
print('minval = ',minval,', maxval = ',maxval)

cm = plt.cm.get_cmap('rainbow')
sc = plt.scatter(distX_SIIa_plot,
                 distY_SIIa_plot,
                 c=np.array(vradSIIa),
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] WIYN')
sb = plt.scatter(np.array(distXaGTC),
                 np.array(distYaGTC),
                 c=np.array(vradaGTC),
                 vmin=minval,
                 vmax=maxval,
                 s=15,
                 cmap=cm,
                 label='[SII] & [ArIII] GTC')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_WIYN_and_GTCa.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

sc = plt.scatter(distX_SIIb_plot,
                 distY_SIIb_plot,
                 c=np.array(vradSIIb),
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] WIYN')
sb = plt.scatter(np.array(distXbGTC),
                 np.array(distYbGTC),
                 c=np.array(vradbGTC),
                 vmin=minval,
                 vmax=maxval,
                 s=15,
                 cmap=cm,
                 label='[SII] GTC')
plt.colorbar(sc, label='radial velocity [km/s]')

plt.ylim = [-110., 110.]
plt.xlim = [-110., 110.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_WIYN_and_GTCb.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()






cm = plt.cm.get_cmap('rainbow')
fig , (ax1,ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12.5,5.0), gridspec_kw = {'width_ratios':[1, 1.27]})
sc = ax1.scatter(distX_SIIa_plot,
                 distY_SIIa_plot,
                 c=np.array(vradSIIa),
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] WIYN')
sb = ax1.scatter(np.array(distXaGTC),
                 np.array(distYaGTC),
                 c=np.array(vradaGTC),
                 vmin=minval,
                 vmax=maxval,
                 s=15,
                 cmap=cm,
                 label='[SII] & [ArIII] GTC')
#plt.colorbar(sc, label='radial velocity [km/s]')

ax1.set_ylim([-110., 110.])
ax1.set_xlim([-110., 110.])
#plt.xlim(-50., -15.)
#plt.xlabel('center distance [arc sec]')
ax1.set_ylabel('center distance [arc sec]')
#ax1.legend()

scb = ax2.scatter(distX_SIIb_plot,
                 distY_SIIb_plot,
                 c=np.array(vradSIIb),
                 vmin=minval,
                 vmax=maxval,
                 s=40,
                 cmap=cm,
                 label='[SII] WIYN')
sbb = ax2.scatter(np.array(distXbGTC),
                 np.array(distYbGTC),
                 c=np.array(vradbGTC),
                 vmin=minval,
                 vmax=maxval,
                 s=15,
                 cmap=cm,
                 label='[SII] GTC')
plt.colorbar(scb, label='radial velocity [km/s]')

ax2.set_ylim([-110., 110.])
ax2.set_xlim([-110., 110.])
#plt.xlim(-50., -15.)
#ax2.set_xlabel('center distance [arc sec]')
#plt.ylabel('center distance [arc sec]')
ax2.legend()
fig.text(0.5,0.02,'center distance [arc sec]',ha='center')
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_vrad_map_WIYN_and_GTC_a+b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()




#sum up spectral range 6690 - 6760 and reproduce an image
wavelengthRange = [6690,6760]
inFile = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_all_fibres_sum.fits'#_6684-6760n.fits'#-skyMean.fits'
#inFile = '/Volumes/work/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMedian.fits'#_6684-6760n.fits'#-skyMean.fits'
hdulist = pyfits.open(inFile)

header = hdulist[0].header
spectrum = pyfits.getdata(inFile)
print('spectrum.shape = ',spectrum.shape)
wavelength = getWavelength(header,1)
print('wavelength.shape = ',wavelength.shape)
print('wavelength = ',wavelength)

SIISum = np.zeros(spectrum.shape[0])
for i in range(spectrum.shape[0]):
    x0 = 0.0
    x1 = 0.0
    for l in range(wavelength.shape[0]):
        if (wavelength[l] >= wavelengthRange[0]) and (wavelength[l] <= wavelengthRange[1]):
            SIISum[i] += spectrum[i,l]
            if x0 == 0:
                x0 = l
            x1 = l
    print('SIISum[',i,'] = ',SIISum[i])
    print('spectrum[',i,',',x0,':',x1,'] = ',spectrum[i,x0:x1])
    print('np.median = ',np.median(spectrum[i,x0:x1]))
    SIISum[i] = SIISum[i] - (np.median(spectrum[i,x0:x1]) * (x1-x0))
    print('SIISum[',i,'] = ',SIISum[i])


print('SIISum = ',SIISum)
SIISumSorted = np.sort(SIISum)
print('SIISumSorted = ',SIISumSorted)
colourRange = [SIISumSorted[6], SIISumSorted[len(SIISumSorted)-6]]

cm = plt.cm.get_cmap('gnuplot')#autumn')#plasma')#copper')
sc = plt.scatter(distX,
                 distY,
                 c=SIISum,
                 vmin=0.,#colourRange[0],
                 vmax=colourRange[1],
                 s=60,
                 cmap=cm)
plt.colorbar(sc, label='ADUs')
plt.title('reconstructed image '+str(wavelengthRange[0])+'-'+str(wavelengthRange[1])+' Ang')

plt.ylim = [-80., 80.]
plt.xlim = [-80., 80.]
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/pa30_reconstructed_image_WIYN.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

from astropy.io import ascii
import numpy as np
import csvFree,csvData
from matplotlib import pyplot as plt


def readCDSTable(tableName,ReadMeName):
    table = ascii.read(tableName,
            readme=ReadMeName)
    return table

clusters = readCDSTable('/Users/azuri/daten/uni/HKU/publications/clusters/clusters.dat','/Users/azuri/daten/uni/HKU/publications/clusters/ReadMe.txt')
if False:
    for cluster in clusters:
        if cluster['Name'] in ['HSC_2686','Lynga_3']:
            print('found cluster. name = ',cluster['Name'])
            print('type = ',cluster['Type']) #o: open cluster, m: moving group, g: globular cluster, d: too distant to classify, r: rejected
            print('Note = ',cluster['Note'])
            print('cluster[N] = ',cluster['N'])
            print('cluster[CST] = ',cluster['CST'])
            print('cluster[RAdeg] = ',cluster['RAdeg'])
            print('cluster[DEdeg] = ',cluster['DEdeg'])
            print('cluster[GLON] = ',cluster['GLON'])
            print('cluster[GLAT] = ',cluster['GLAT'])
            print('cluster[pmRA] = ',cluster['pmRA'],' +/- ',cluster['e_pmRA'])
            print('cluster[pmDE] = ',cluster['pmDE'],' +/- ',cluster['e_pmDE'])
            print('cluster[s_pmRA] = ',cluster['s_pmRA'])
            print('cluster[s_pmDE] = ',cluster['s_pmDE'])
            print('cluster[Plx] = ',cluster['Plx'],' +/- ',cluster['e_Plx'])
            print('cluster[s_Plx] = ',cluster['s_Plx'])
            print('cluster[dist16] = ',cluster['dist16'])
            print('cluster[dist50] = ',cluster['dist50'])
            print('cluster[dist84] = ',cluster['dist84'])
            print('cluster[X] = ',cluster['X'])
            print('cluster[Y] = ',cluster['Y'])
            print('cluster[Z] = ',cluster['Z'])
            print('cluster[RV] = ',cluster['RV'])
            print('cluster[e_RV] = ',cluster['e_RV'])
            print('cluster[s_RV] = ',cluster['s_RV'])
            print('cluster[n_RV] = ',cluster['n_RV'])
            print('cluster[logAge16] = ',cluster['logAge16'])
            print('cluster[logAge50] = ',cluster['logAge50'])
            print('cluster[logAge84] = ',cluster['logAge84'])
            print('cluster[AV16] = ',cluster['AV16'])
            print('cluster[AV50] = ',cluster['AV50'])
            print('cluster[AV84] = ',cluster['AV84'])
            print('cluster[CMDCl2.5] = ',cluster['CMDCl2.5'])
            print('cluster[CMDCl16] = ',cluster['CMDCl16'])
            print('cluster[CMDCl50] = ',cluster['CMDCl50'])
            print('cluster[CMDCl84] = ',cluster['CMDCl84'])
            print('cluster[CMDCl97.5] = ',cluster['CMDCl97.5'])
            print('cluster[CMDClHuman] = ',cluster['CMDClHuman'])

members = readCDSTable('/Users/azuri/daten/uni/HKU/publications/clusters/members.dat','/Users/azuri/daten/uni/HKU/publications/clusters/ReadMe.txt')
print('dir(members) = ',dir(members))
print('members.colnames = ',members.colnames)
#STOP
if True:
    hsc_RA = []
    hsc_DEC = []
    hsc_pmRA = []
    hsc_pmDEC = []
    hsc_epmRA = []
    hsc_epmDEC = []

    lynga_RA = []
    lynga_DEC = []
    lynga_pmRA = []
    lynga_pmDEC = []
    lynga_epmRA = []
    lynga_epmDEC = []

    cs_RA = []
    cs_DEC = []
    cs_pmRA = []
    cs_pmDEC = []
    cs_epmRA = []
    cs_epmDEC = []

    for member in members:
        #print('type(member[GaiaDR3]) = ',type(member['GaiaDR3']))
        #STOP
        if member['Name'] == 'HSC_2686':
            hsc_RA.append(member['RAdeg'])
            hsc_DEC.append(member['DEdeg'])
            hsc_pmRA.append(member['pmRA'])
            hsc_pmDEC.append(member['pmDE'])
            hsc_epmRA.append(member['e_pmRA'])
            hsc_epmDEC.append(member['e_pmDE'])
        elif member['Name'] == 'Lynga_3':
            lynga_RA.append(member['RAdeg'])
            lynga_DEC.append(member['DEdeg'])
            lynga_pmRA.append(member['pmRA'])
            lynga_pmDEC.append(member['pmDE'])
            lynga_epmRA.append(member['e_pmRA'])
            lynga_epmDEC.append(member['e_pmDE'])

        if member['GaiaDR3'] == 5877088151005347072:
            print('member[Name] = ',member['Name'])
            print('member[RAdeg] = ',member['RAdeg'])
            print('member[DEdeg] = ',member['DEdeg'])
            print('member[pmRA] = ',member['pmRA'],'+/-',member['e_pmRA'])
            print('member[pmDE] = ',member['pmDE'],'+/-',member['e_pmDE'])
            print('member[Plx] = ',member['Plx'],'+/-',member['e_Plx'])
            print('member[RV] = ',member['RV'],'+/-',member['e_RV'])
            cs_RA.append(member['RAdeg'])
            cs_DEC.append(member['DEdeg'])
            cs_pmRA.append(member['pmRA'])
            cs_pmDEC.append(member['pmDE'])
            cs_epmRA.append(member['e_pmRA'])
            cs_epmDEC.append(member['e_pmDE'])

gaia = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/julia/gaia_all_stars_in_area.csv')
print('gaia.header = ',gaia.header)

sourceID = np.array(gaia.getData('SOURCE_ID'))
ra = np.asarray(gaia.getData('ra'),dtype=float)
dec = np.asarray(gaia.getData('dec'),dtype=float)
print('ra = ',ra)
print('dec = ',dec)
plt.scatter(ra,dec,s=0.5,c='grey')
plt.scatter(hsc_RA,hsc_DEC,s=2,c='r',label='HSC 2686')
plt.scatter(lynga_RA,lynga_DEC,s=2,c='b',label='Lynga 3')

pmra = np.array(gaia.getData('pmra'))
print('pmra = ',pmra)
idx = np.where(pmra != '')[0]
print('idx = ',idx)
pmra = pmra[idx]
pmdec = np.array(gaia.getData('pmdec'))[idx]
pmra = np.asarray(pmra,dtype=float)
pmdec = np.asarray(pmdec,dtype=float)

idx = np.where(pmra < -5.2)[0]
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_ID = sourceID[idx]

idx = np.where(rand_pmra > -5.4)
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_ID = sourceID[idx]

idx = np.where(rand_pmdec < -3.8)
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_ID = sourceID[idx]

idx = np.where(rand_pmdec > -4.)
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_ID = sourceID[idx]

idx = np.random.randint(0,len(rand_ID),len(hsc_DEC))
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_ID = sourceID[idx]

plt.scatter(rand_ra,rand_dec,s=5,c='m',label='random sample')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.legend()
plt.show()



plt.scatter(pmra,pmdec,s=0.5,c='grey')
plt.scatter(hsc_pmRA,hsc_pmDEC,s=2,c='r',label='HSC 2686')
plt.errorbar(hsc_pmRA,hsc_pmDEC, xerr=hsc_epmRA, yerr=hsc_epmDEC, fmt="o", color="r")
plt.scatter(lynga_pmRA,lynga_pmDEC,s=2,c='b',label='Lynga 3')
plt.errorbar(lynga_pmRA,lynga_pmDEC, xerr=lynga_epmRA, yerr=lynga_epmDEC, fmt="o", color="b")
plt.scatter(cs_pmRA,cs_pmDEC,s=5,c='g',label='CSPN')
plt.errorbar(cs_pmRA,cs_pmDEC, xerr=cs_epmRA, yerr=cs_epmDEC, fmt="o", color="g")
plt.scatter(rand_pmra,rand_pmdec,s=10,c='m',label='random sample')
plt.xlim(-5.6,-4.5)
plt.ylim(-4.4,-3)
plt.xlabel('pm_RA')
plt.ylabel('pm_DEC')
plt.legend()
plt.show()



parallax = np.array(gaia.getData('parallax'))
idx = np.where(parallax != '')[0]
parallax = parallax[idx]
vrad = np.array(gaia.getData('radial_velocity'))[idx]
idx = np.where(vrad != '')[0]
parallax = np.asarray(parallax[idx],dtype=float)
vrad = np.asarray(vrad[idx],dtype=float)
plt.scatter(vrad,parallax,s=0.5)
plt.xlabel('radial velocity')
plt.ylabel('parallax')
plt.show()

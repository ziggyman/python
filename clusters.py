from astropy.io import ascii
import numpy as np
import csvFree,csvData
from matplotlib import pyplot as plt


def readCDSTable(tableName,ReadMeName):
    table = ascii.read(tableName,
            readme=ReadMeName)
    return table


def distance(plx):
    return 1000. / plx  # mas to pc

def apparant_Gmag(Gmag, plx):
    return 5 * np.log10(distance(plx)) - 5 + Gmag

def createGaiaIsoTable():
    import glob
    import os
    from numpy.lib import recfunctions as rfn
    files = [f for f in sorted(glob.glob("/Users/azuri/daten/uni/HKU/isochrones/*.dat", recursive=True),key=os.path.getsize)]
    print('files = ',files)
    sarr = None
    for file in files:
        if sarr is None:
            sarr = np.genfromtxt(file,comments='#',skip_header=14,names=True)
        else:
            sarr = rfn.stack_arrays((sarr,np.genfromtxt(file,comments='#',skip_header=14,names=True)))
        print('file = ',file,': ',sarr.shape,': sarr.names = ',sarr.dtype.names)
    #np.save('/Users/azuri/entwicklung/python/OCFit/gaiaDR2/grids/full_isoc_GAIA_CMD33.npy',sarr, allow_pickle=True)
    #sarr = np.load('/Users/azuri/entwicklung/python/OCFit/gaiaDR2/grids/full_isoc_GAIA_CMD33.npy')
    #print('sarr = ',sarr)
    return sarr

#createGaiaIsoTable()
#STOP

def createNewCSVDataFromIndices(csvIn,idxIn):
    csvOut = csvData.CSVData()
    csvOut.header = csvIn.header
    for idx in idxIn:
        csvOut.append(csvIn.getData(idx))
    return csvOut

clusters = readCDSTable('/Users/azuri/daten/uni/HKU/publications/clusters/clusters.dat','/Users/azuri/daten/uni/HKU/publications/clusters/ReadMe.txt')
if True:
    for cluster in clusters:
        if cluster['Name'] in ['NGC_6530','Trumpler_14','NGC_2244','HSC_2686','Lynga_3',]:
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
#STOP
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
    hsc_Gmag = []
    hsc_BPmag = []
    hsc_RPmag = []
    hsc_parallax = []
    hsc_e_parallax = []

    lynga_RA = []
    lynga_DEC = []
    lynga_pmRA = []
    lynga_pmDEC = []
    lynga_epmRA = []
    lynga_epmDEC = []
    lynga_Gmag = []
    lynga_BPmag = []
    lynga_RPmag = []
    lynga_parallax = []
    lynga_e_parallax = []

    cs_RA = []
    cs_DEC = []
    cs_pmRA = []
    cs_pmDEC = []
    cs_epmRA = []
    cs_epmDEC = []
    cs_Gmag = []
    cs_BPmag = []
    cs_RPmag = []
    cs_parallax = []
    cs_e_parallax = []

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
            hsc_Gmag.append(member['Gmag'])
            hsc_BPmag.append(member['BPmag'])
            hsc_RPmag.append(member['RPmag'])
            hsc_parallax.append(member['Plx'])
            hsc_e_parallax.append(member['e_Plx'])
        elif member['Name'] == 'Lynga_3':
            lynga_RA.append(member['RAdeg'])
            lynga_DEC.append(member['DEdeg'])
            lynga_pmRA.append(member['pmRA'])
            lynga_pmDEC.append(member['pmDE'])
            lynga_epmRA.append(member['e_pmRA'])
            lynga_epmDEC.append(member['e_pmDE'])
            lynga_Gmag.append(member['Gmag'])
            lynga_BPmag.append(member['BPmag'])
            lynga_RPmag.append(member['RPmag'])
            lynga_parallax.append(member['Plx'])
            lynga_e_parallax.append(member['e_Plx'])

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
            cs_Gmag.append(member['Gmag'])
            cs_BPmag.append(member['BPmag'])
            cs_RPmag.append(member['RPmag'])
            cs_parallax.append(member['Plx'])
            cs_e_parallax.append(member['e_Plx'])

    hsc_Gmag = np.asarray(hsc_Gmag,dtype=float)
    hsc_BPmag = np.asarray(hsc_BPmag,dtype=float)
    hsc_RPmag = np.asarray(hsc_RPmag,dtype=float)
    hsc_parallax = np.asarray(hsc_parallax,dtype=float)
    hsc_e_parallax = np.asarray(hsc_e_parallax,dtype=float)
    hsc_pmRA = np.asarray(hsc_pmRA,dtype=float)
    hsc_epmRA = np.asarray(hsc_epmRA,dtype=float)
    hsc_pmDEC = np.asarray(hsc_pmDEC,dtype=float)
    hsc_epmDEC = np.asarray(hsc_epmDEC,dtype=float)

    lynga_Gmag = np.asarray(lynga_Gmag,dtype=float)
    lynga_BPmag = np.asarray(lynga_BPmag,dtype=float)
    lynga_RPmag = np.asarray(lynga_RPmag,dtype=float)
    lynga_parallax = np.asarray(lynga_parallax,dtype=float)
    lynga_e_parallax = np.asarray(lynga_e_parallax,dtype=float)
    lynga_pmRA = np.asarray(lynga_pmRA,dtype=float)
    lynga_epmRA = np.asarray(lynga_epmRA,dtype=float)
    lynga_pmDEC = np.asarray(lynga_pmDEC,dtype=float)
    lynga_epmDEC = np.asarray(lynga_epmDEC,dtype=float)

    cs_Gmag = np.asarray(cs_Gmag,dtype=float)
    cs_BPmag = np.asarray(cs_BPmag,dtype=float)
    cs_RPmag = np.asarray(cs_RPmag,dtype=float)
    cs_parallax = np.asarray(cs_parallax,dtype=float)
    cs_e_parallax = np.asarray(cs_e_parallax,dtype=float)
    cs_pmRA = np.asarray(cs_pmRA,dtype=float)
    cs_epmRA = np.asarray(cs_epmRA,dtype=float)
    cs_pmDEC = np.asarray(cs_pmDEC,dtype=float)
    cs_epmDEC = np.asarray(cs_epmDEC,dtype=float)

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
pmdec = np.array(gaia.getData('pmdec'))
parallax = np.array(gaia.getData('parallax'))
Gmag = np.array(gaia.getData('phot_g_mean_mag'))
BPmag = np.array(gaia.getData('phot_bp_mean_mag'))
RPmag = np.array(gaia.getData('phot_rp_mean_mag'))
print('pmra = ',pmra)
cond1 = pmra != ''
cond2 = parallax != ''
cond3 = Gmag != ''
cond4 = BPmag != ''
cond5 = RPmag != ''
idx = np.asarray(np.where(cond1 & cond2 & cond3 & cond4 & cond5)[0],dtype=int)
print('idx = ',idx)
pmra = np.asarray(pmra[idx],dtype=float)
pmdec = np.asarray(pmdec[idx],dtype=float)
sourceID = sourceID[idx]
ra = ra[idx]
dec = dec[idx]
Gmag = Gmag[idx]
Gmag = np.asarray(Gmag,dtype=float)
#Gmag = np.asarray(Gmag[idx],dtype=float)
BPmag = np.asarray(BPmag[idx],dtype=float)
print('idx = ',type(idx),': ',type(idx[0]),': ',idx)
RPmag = np.asarray(RPmag[idx],dtype=float)
parallax = np.asarray(parallax[idx],dtype=float)
gaiaOut = createNewCSVDataFromIndices(gaia,idx)

print('parallax.shape = ',parallax.shape,', pmra.shape = ',pmra.shape,', pmdec.shape = ',pmdec.shape)

cond1 = pmra < -5.2
cond2 = pmra > -5.4
cond3 = pmdec < -3.8
cond4 = pmdec > -4.
cond5 = parallax > 0.08
cond6 = parallax < 0.2

idx = np.where(cond1 & cond2 & cond3 & cond4 & cond5 & cond6)[0]
print('random: len(idx) = ',len(idx))

rand_ID = sourceID[idx]
rand_ra = ra[idx]
rand_dec = dec[idx]
rand_pmra = pmra[idx]
rand_pmdec = pmdec[idx]
rand_Gmag = Gmag[idx]
rand_BPmag = BPmag[idx]
rand_RPmag = RPmag[idx]
rand_parallax = parallax[idx]
gaiaOut = createNewCSVDataFromIndices(gaiaOut,idx)

if False:
    idx = np.random.randint(0,len(rand_ID),len(hsc_DEC))
    print('idx = ',len(idx),': ',idx)
    rand_ra = rand_ra[idx]
    rand_dec = rand_dec[idx]
    rand_pmra = rand_pmra[idx]
    rand_pmdec = rand_pmdec[idx]
    rand_ID = rand_ID[idx]
    rand_Gmag = rand_Gmag[idx]
    rand_BPmag = rand_BPmag[idx]
    rand_RPmag = rand_RPmag[idx]
    rand_parallax = rand_parallax[idx]
    gaiaOut = createNewCSVDataFromIndices(gaiaOut,idx)
print('rand_pmra = ',rand_pmra)
print('rand_pmdec = ',rand_pmdec)
#STOP
gaiaOut.renameColumn('ra','RAdeg')
gaiaOut.renameColumn('ra_error','e_RAdeg')
gaiaOut.renameColumn('dec','DEdeg')
gaiaOut.renameColumn('dec_error','e_DEdeg')
gaiaOut.renameColumn('parallax','Plx')
gaiaOut.renameColumn('parallax_error','e_Plx')
gaiaOut.renameColumn('pmra','pmRA')
gaiaOut.renameColumn('pmra_error','e_pmRA')
gaiaOut.renameColumn('pmdec','pmDE')
gaiaOut.renameColumn('pmdec_error','e_pmDE')
gaiaOut.renameColumn('phot_g_mean_flux','FG')
gaiaOut.renameColumn('phot_g_mean_flux_error','e_FG')
gaiaOut.renameColumn('phot_bp_mean_flux','FBP')
gaiaOut.renameColumn('phot_bp_mean_flux_error','e_FBP')
gaiaOut.renameColumn('phot_rp_mean_flux','FRP')
gaiaOut.renameColumn('phot_rp_mean_flux_error','e_FRP')
gaiaOut.renameColumn('phot_g_mean_mag','Gmag')
#gaiaOut.renameColumn('phot_g_mean_mag_error','e_Gmag')
gaiaOut.renameColumn('phot_bp_mean_mag','BPmag')
gaiaOut.renameColumn('phot_rp_mean_mag','RPmag')
gaiaOut.renameColumn('radial_velocity','RV')
gaiaOut.renameColumn('radial_velocity_error','e_RV')
gaiaOut.addColumn('Prob',np.zeros(gaiaOut.size())+0.9)
gaiaOut.addColumn('Name',['random_sample3' for i in range(gaiaOut.size())])
csvFree.writeCSVFile(gaiaOut,'/Users/azuri/daten/uni/HKU/publications/clusters/gaia_random_sample3.csv')

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
#plt.xlim(-5.6,-4.5)
#plt.ylim(-4.4,-3)
plt.xlabel('pm_RA')
plt.ylabel('pm_DEC')
plt.legend()
plt.show()


plt.scatter(np.sqrt(pmra**2 + pmdec**2), parallax,s=0.5,c='grey')
plt.scatter(np.sqrt(hsc_pmRA**2 + hsc_pmDEC**2), hsc_parallax,s=2,c='r',label='HSC 2686')
plt.scatter(np.sqrt(lynga_pmRA**2 + lynga_pmDEC**2), lynga_parallax,s=2,c='b',label='Lynga 3')
plt.scatter(np.sqrt(rand_pmra**2+rand_pmdec**2),rand_parallax,s=10,c='m',label='random sample')
plt.xlabel('pm [mas/yr]')
plt.ylabel('parallax [mas]')
plt.legend()
plt.show()



vrad = np.array(gaia.getData('radial_velocity'))[idx]
idx = np.where(vrad != '')[0]
parallax = np.asarray(parallax[idx],dtype=float)
vrad = np.asarray(vrad[idx],dtype=float)
plt.scatter(vrad,parallax,s=0.5)
plt.xlabel('radial velocity')
plt.ylabel('parallax')
plt.show()



plt.scatter(hsc_BPmag-hsc_RPmag,hsc_Gmag,c='r',s=5,label='HSC 2686')
plt.scatter(lynga_BPmag-lynga_RPmag,lynga_Gmag,c='b',s=5,label='Lynga 3')
plt.scatter(rand_BPmag-rand_RPmag,rand_Gmag,c='m',s=5,label = 'random sample')
plt.xlabel('G_BP-G_RP [mag]')
plt.ylabel('G [mag]')
plt.legend()
plt.show()

iso = np.genfromtxt('/Users/azuri/daten/uni/HKU/isochrones/HSC_2686_from_HR_best_fit.dat',comments='#',skip_header=14,names=True)
metallicity_grid = np.unique(iso['MH'])
print('metallicity_grid = ',metallicity_grid)
age_grid = np.unique(iso['logAge'])
print('age_grid = ',age_grid)
dist = distance(0.15516932)
print('dist = ',dist)

plt.scatter(hsc_BPmag-hsc_RPmag,hsc_Gmag,c='r',s=5,label='HSC 2686')
for iMH in range(len(metallicity_grid)):
    isochrone = iso[np.where(iso['MH'] == metallicity_grid[iMH])]
    #plt.scatter(apparant_Gmag(isochrone['G_BPmag'],dist) - apparant_Gmag(isochrone['G_RPmag'],dist), apparant_Gmag(isochrone['Gmag'],dist),s=2)
    plt.scatter(isochrone['G_BPmag'] - isochrone['G_RPmag'], isochrone['Gmag'],s=2)
plt.xlabel('G_BP-G_RP [mag]')
plt.ylabel('G [mag]')
plt.legend()
plt.show()

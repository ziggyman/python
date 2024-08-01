from drUtils import traceMultiApertureImage,fitEmissionLines

spec='/Users/azuri/spectra/MSO/MSSO_Feb2007/MEN0844-4600_MS220207.fits'
fitEmissionLines(spec,'/Users/azuri/Downloads/MEN0844-4600_MS220207.fits')
STOP
traceMultiApertureImage('/Users/azuri/spectra/saao/saao_nov-2016/saao_nov2016/data/refVerticalTrace_spupnic_gr7_15_7_Nov2013_otzf.fits','/Users/azuri/spectra/saao/saao_nov-2016/saao_nov2016/data/database/aprefVerticalTrace_spupnic_gr7_15_7_Nov2013_otzf.new',peakHeight=400,peakWidth=2.,threshold=100,step=1,nSum=10,nLost=3,nPeaksShould=None,transpose=True)

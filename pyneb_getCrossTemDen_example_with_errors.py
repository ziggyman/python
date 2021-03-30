import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

# Tell PyNeb tu use parallelisation
pn.config.use_multiprocs()        

### General settings
# Setting verbosity level. Enter pn.my_logging? for details
pn.log_.level = 2 # set this to 3 to have more details

# Adopt an extinction law
extinction_law = 'CCM89'

# Define the data file
obs_data = 'test.dat'

# Define plot title
title = 'SMC 24'

### Read and deredden observational data
# define an Observation object and assign it to name 'obs'
obs = pn.Observation()

# fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
obs.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)

# Add a number of "fake" observations aroung the right one, to make Monte Carlo statistics
obs.addMonteCarloObs(N = 500)

obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
# deredden data with Cardelli's law
obs.extinction.law = extinction_law
obs.correctData()

print("PyNeb: cH(b) = {:.2f} +/- {:.2f}".format(np.median(obs.extinction.cHbeta),np.std(obs.extinction.cHbeta)))
### Include the diagnostics of interest
# instantiate the Diagnostics class
diags = pn.Diagnostics()
# include in diags the relevant line ratios
diags.addDiag([
              '[NII] 5755/6584', 
              '[OII] 3726/3729', 
              '[OIII] 4363/5007', 
              '[SII] 6731/6716', 
              '[SII] 4072+/6720+',
              '[SIII] 6312/18.7m', 
              '[NeIII] 3930+/15.6m', 
              ])
diags.addClabel('[SII] 6731/6716', '[SII]a')
diags.addClabel('[SII] 4072+/6720+', '[SII]b')



#The observed ratio can be automatically extracted from an Observation object named obs:
Te, Ne = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716', obs=obs)
print('Te = {:.0f} +/- {:.0f}'.format(np.mean(Te), np.std(Te)))
print('Ne = {:.0f} +/- {:.0f}'.format(np.mean(Ne), np.std(Ne)))



# creating a new .dat file to let the pyneb lib to read the data
# file structure: for each object, there is a file, in each file, each row represent 1 row in the spectrum

for i in range(len(object_name)):
    with open('LSR//intensity//updated//observation_' + str(i) + '.dat', mode='w') as f:
        f.writelines('H1r_4861A\tH1r_4861Ae\tH1r_6563A\tH1r_6563Ae\tN2_5755A\tN2_5755Ae\tN2_6548A\tN2_6548Ae\tN2_6584A\tN2_6584Ae\tO3_4363A\tO3_4363Ae\tO3_5007A\tO3_5007Ae\tS2_6716A\tS2_6716Ae\tS2_6731A\tS2_6731Ae\tAr4_4740A\tAr4_4740Ae\tAr4_4711A\tAr4_4711Ae\tAr5_4626A\tAr5_4626Ae\tAr5_6435A\tAr5_6435Ae\tAr5_7005A\tAr5_7005Ae\n')
        for j in range(90):
            f.writelines(str(intensity[i,j,10])+'\t'+str(intensity[i,j,10] * 0.1)+'\t'+str(intensity[i,j,34])+'\t'+str(intensity[i,j,34] * 0.1)+'\t'
                    +str(intensity[i,j,26])+'\t'+str(intensity[i,j,26] * 0.1)+'\t'+str(intensity[i,j,33])+'\t'+str(intensity[i,j,33] * 0.1)+'\t'
                    +str(intensity[i,j,36])+'\t'+str(intensity[i,j,36] * 0.1)+'\t'+str(intensity[i,j,63])+'\t'+str(intensity[i,j,63] * 0.1)+'\t'
                    +str(intensity[i,j,62])+'\t'+str(intensity[i,j,62] * 0.1)+'\t'+str(intensity[i,j,38])+'\t'+str(intensity[i,j,38] * 0.1)+'\t'
                    +str(intensity[i,j,39])+'\t'+str(intensity[i,j,39] * 0.1)+'\t'+str(intensity[i,j,5])+'\t'+str(intensity[i,j,5] * 0.1)+'\t'
                    +str(intensity[i,j,8])+'\t'+str(intensity[i,j,8] * 0.1)+'\t'+str(intensity[i,j,64])+'\t'+str(intensity[i,j,64] * 0.1)+'\t'
                    +str(intensity[i,j,65])+'\t'+str(intensity[i,j,65] * 0.1)+'\t'+str(intensity[i,j,66])+'\t'+str(intensity[i,j,66] * 0.1)+'\n')
    f.close()



import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

extinction = np.empty([len(object_name),90])
ext_err = np.empty(len(object_name))
extinction_all = []

fig_ext = []
fig_temp_OIII = []
fig_density_OIII = []
fig_temp_ArV = []
fig_density_ArV = []
temp_OIII_err = []
density_OIII_err = []

# Tell PyNeb tu use parallelisation
pn.config.use_multiprocs()

### General settings
# Setting verbosity level. Enter pn.my_logging? for details
pn.log_.level = 2 # set this to 3 to have more details

# Adopt an extinction law
extinction_law = 'CCM89'

index = [0,4,6,8,9,10,11,12,14,20,25]

for j in range (len(index)):
# i = 6
    i = index[j]
    ### Read and deredden observational data
    # define an Observation object and assign it to name 'obs'
    obs = pn.Observation()

    #fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData('LSR//intensity//updated//observation_' + str(i) + '.dat', errIsRelative=False)# fileFormat='lines_in_rows_err_cols', err_default=0.05)

    obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
    # deredden data with Cardelli's law
    obs.extinction.law = extinction_law
    obs.correctData()

    extinction[i] = obs.extinction.cHbeta

    # Plot the filtered data
    temp_cen = cen[i]
    temp_ini = (temp_cen + 1) * (-0.9)
    temp_fin = (69 - temp_cen) * 0.9
    x_axis = np.linspace(temp_ini, temp_fin, 90)

    fig, ax = plt.subplots()
    ax.scatter(x_axis, extinction[i])
    ax.set_ylabel("Extinction")
    ax.set_xlabel("ArcSec")
    ax.set_xlim(temp_ini, temp_fin)
    ax.set_title("Object: " + str(object_name[i]))
    fig_ext.append(fig)


    # finding the temperature and density
    ### Include the diagnostics of interest
    # instantiate the Diagnostics class
    diags = pn.Diagnostics()

    #print(obs.printIntens())


    # include in diags the relevant line ratios
    diags.addDiag([
                  '[NII] 5755/6584',
                  '[OII] 3726/3729',
                  '[OIII] 4363/5007',
                  '[SII] 6731/6716',
                  '[ArIV] 4740/4711',
                  '[ArV] 4626/6600+',
                  '[SII] 4072+/6720+',
                  '[SIII] 6312/18.7m',
                  '[NeIII] 3930+/15.6m',
                  ])
    diags.addClabel('[SII] 6731/6716', '[SII]a')
    diags.addClabel('[SII] 4072+/6720+', '[SII]b')



    #The observed ratio can be automatically extracted from an Observation object named obs:
    print(i)
    print("Using OIII for temperature")
    Te_O, Ne_O = diags.getCrossTemDen(diag_tem='[OIII] 4363/5007', diag_den='[SII] 6731/6716', obs = obs)
                                  #value_tem=(o_4363[i]/o_5007[i]), value_den=(s_6731[i]/s_6716[i]))
    print("temp = ", Te_O)
    print("den = ", Ne_O)



    # plotting the graphs
    temp_cen = cen[i]
    temp_ini = (temp_cen + 1) * (-0.9)
    temp_fin = (69 - temp_cen) * 0.9
    x_axis = np.linspace(temp_ini, temp_fin, 90)


    fig1, ax1 = plt.subplots()
    ax1.scatter(x_axis, Ne_O, label="Electron Density", color='orange')
    ax1.set_ylabel("Electron Density (cm^-3)")
    ax1.set_xlabel("Arcsec")
    #ax1.set_xlim(temp_ini, temp_fin)
    ax1.set_xlim(-25,25)
    #ax2.set_ylim(0, np.nanmax(temp_electron_density) * 1.1)
    ax1.set_title(f"Electron Density for Object: {object_name[i]} Using OIII")
    ax1.legend()
    fig_density_OIII.append(fig1)

    fig3, ax3 = plt.subplots()
    ax3.scatter(x_axis, Te_O, label="Temperature", color='orange')
    ax3.set_ylabel("Temperature (K)")
    ax3.set_xlabel("Arcsec")
    #ax3.set_xlim(temp_ini, temp_fin)
    ax3.set_xlim(-25,25)
    #ax4.set_ylim(0, np.nanmax(temp_temperature) * 1.1)
    ax3.set_title(f"Temperature for Object: {object_name[i]} Using OIII")
    ax3.legend()
    fig_temp_OIII.append(fig3)


    # calculating the error

    obs = pn.Observation()

    #fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData('LSR//intensity//observation_' + str(i) + '.dat')# fileFormat='lines_in_rows_err_cols', errIsRelative=False, err_default=0.05)

    # Add a number of "fake" observations aroung the right one, to make Monte Carlo statistics
    obs.addMonteCarloObs(N = 500) # finding the error which is the sd of the data set

    obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
    # deredden data with Cardelli's law
    obs.extinction.law = extinction_law
    obs.correctData()

    # create a new set of data from obs.extinction.cHbeta which only contain useful data but not nan/inf
    temp_ext = []
    for m in range(len(obs.extinction.cHbeta)):
        if obs.extinction.cHbeta[m] > 0 and obs.extinction.cHbeta[m] < 10000:
            temp_ext.append(obs.extinction.cHbeta[m])

    ext_err[i] = np.std(temp_ext)

    print("extinction error = ", ext_err[i])



    #The observed ratio can be automatically extracted from an Observation object named obs:
    print("calculating the error")
    Te, Ne = diags.getCrossTemDen(diag_tem='[OIII] 4363/5007', diag_den='[SII] 6731/6716', obs = obs)
                                  #value_tem=(o_4363[i]/o_5007[i]), value_den=(s_6731[i]/s_6716[i]))

    temp_temp_err = []
    temp_density_err = []

    for l in range (90):
        num = l*500
        temp_temp_err.append(np.std(Te[(0+num):(499+num)]))
        temp_density_err.append(np.std(Ne[(0+num):(499+num)]))


    temp_OIII_err.append(temp_temp_err)
    density_OIII_err.append(temp_density_err)

    print("temp err: ", temp_temp_err)
    print("density err: ", temp_density_err)

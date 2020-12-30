import subprocess

# Input data
# Specify either l,b or rlen. Set other to NA. rlen takes precedence.
w    = 1 # parallax in mas (corrected for any zeropoint offset; +0.029mas in the catalogue)
wsd  = abs(0.2*w) # parallax uncertainty in mas
glon = 340 # Galactic longitude in degrees (0 to 360)
glat =  45 # Galactic latitude (-90 to +90)
rlen =  "NA" # length scale in pc
# Plotting parameters in pc
# rlo,rhi are range for computing normalization of posterior
# rplotlo, rplothi are plotting range (computed automatically if set to NA)
rlo = 0
rhi = 1e5
rplotlo = "NA"
rplothi = "NA"

# Define command and arguments
command ='Rscript'
path2script ='/Users/azuri/entwicklung/r/Gaia-DR2-distances/run_distest_single.R'

# Variable number of args in a list
args = [str(w),str(wsd),str(glon),str(glat),rlen,str(rlo),str(rhi),rplotlo,rplothi]

# Build subprocess command
cmd = [command, path2script] + args

# check_output will run the command and store to result
x = subprocess.check_output(cmd, universal_newlines=True)

print('The maximum of the numbers is:', x)

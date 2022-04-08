import matplotlib.pyplot as plt
import os
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

listName = '/Users/azuri/spectra/saao/saao_may2007/RAW/night5/FLAT.list'
with open(listName,'r') as f:
    lines = f.readlines()
lines = [line.strip() for line in lines]

for fName in lines:
    #fName = os.path.join(listName[:listName.rfind('/')],line)
    print('fName = <'+fName+'>')
    #image_file = get_pkg_data_filename(fName)
    #fits.info(image_file)
    #image_data = fits.getdata(image_file, ext=0)
    hdu = fits.open(fName)[0]
    image_data = hdu.data
    print(image_data.shape)
    plt.figure()
    plt.imshow(image_data, cmap='gray')
    #plt.colorbar()
    plt.savefig(fName[:fName.rfind('.')]+'.png', bbox_inches='tight')
    plt.show()

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

tab = [
#{'fiber': 1135, 'centerDistanceX' : -85.12, 'centerDistanceY' : -85.12, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6724.81, 'SII6731b' : 6737.28, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'fiber': 1, 'centerDistanceX' : 9.85, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6696.62, 18289.,20.9], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6754.51, 10737., 23.6], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 2, 'centerDistanceX' : -29.56, 'centerDistanceY' : 61.92, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6713.3,7149.0,5.672], 'SII6731a' : [6726.07,7619.0,7.832], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 3, 'centerDistanceX' : 34.49, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6730.48,3930.0,5.12], 'SII6731b' : [6745.08,6399.0,6.95], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 4, 'centerDistanceX' : 24.64, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6703.77,4931.0,8.954], 'SII6731a' : [6721.79,9516.0,13.13], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 5, 'centerDistanceX' : 29.56, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6716.82,2505.0,4.993], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6746.48,12718.0,21.9], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 6, 'centerDistanceX' : 29.56, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.26,76862.0,9.071], 'SII6731a' : [6717.98,7442.0,9.852], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 7, 'centerDistanceX' : 24.64, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6693.59,5645.0,13.32], 'SII6731a' : [6711.32,8018.0,15.05], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 8, 'centerDistanceX' : 29.56, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6700.25,4676.0,11.78], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 9, 'centerDistanceX' : 34.49, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 10, 'centerDistanceX' : 34.49, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6726.09,6722.0,11.31], 'SII6731b' : [6743.01,10599.0,14.66], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 11, 'centerDistanceX' : 24.64, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6698.94,6805.0,8.99], 'SII6731a' : [6714.9,3550.0,4.009], 'SII6716b' : [6730.94,4552.0,7.866], 'SII6731b' : [6743.89,8403.0,13.03], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 12, 'centerDistanceX' : 29.56, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6704.1,4780.0,7.665], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6745.27,7437.0,10.23], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 13, 'centerDistanceX' : 34.49, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6727.22,4436.0,10.62], 'SII6731b' : [6747.52,4568.0,6.754], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 14, 'centerDistanceX' : 29.56, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6739.42,6693.0,6.918], 'SII6731b' : [6753.71,4630.0,5.577], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 15, 'centerDistanceX' : 24.64, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6696.13,3748.0,5.126], 'SII6731a' : [6708.96,4156.0,10.11], 'SII6716b' : [6730.28,4694.0,6.399], 'SII6731b' : [6744.81,7185.0,8.379], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 16, 'centerDistanceX' : 0.0, 'centerDistanceY' : 61.92, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6746.04,5866.0,8.922], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 17, 'centerDistanceX' : 19.71, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.74,4785.0,7.129], 'SII6731a' : [6722.02,10062.0,14.37], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6744.92,3591.0,5.948], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 18, 'centerDistanceX' : 4.93, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [4988.85,20479.0,2.921], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6695.67,9157.0,7.79], 'SII6731a' : [6708.06,8492.0,8.629], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 19, 'centerDistanceX' : 9.85, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6701.8,3261.0,3.222], 'SII6731a' : [6714.96,6414.0,10.71], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 20, 'centerDistanceX' : 14.78, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6703.34,4708.0,6.392], 'SII6731a' : [6717.12,6281.0,7.234], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 21, 'centerDistanceX' : 9.85, 'centerDistanceY' : 16.89, 'Halpha6563a': [6538.0,5087.0,2.931], 'Halpha6563b': [6581.12,4449.0,3.129], 'OIII5007a': [4990.98,27483.0,3.864], 'OIII5007b': [5020.9,15711.0,2.96], 'SII6716a': [6696.75,17792.0,6.442], 'SII6731a' : [6711.58,14306.0,7.242], 'SII6716b' : [6734.78,5993.0,7.468], 'SII6731b' : [6751.43,4412.0,8.583], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 22, 'centerDistanceX' : 29.56, 'centerDistanceY' : 61.92, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6713.41,2359.0,3.665], 'SII6731a' : [6725.49,7284.0,7.063], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6744.61,3890.0,6.373], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 23, 'centerDistanceX' : 14.78, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6706.16,4161.0,7.101], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 24, 'centerDistanceX' : 14.78, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 25, 'centerDistanceX' : 4.93, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 26, 'centerDistanceX' : 4.93, 'centerDistanceY' : 2.81, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 27, 'centerDistanceX' : 19.71, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6699.54,5106.0,8.047], 'SII6731a' : [6714.31,2995.0,3.896], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 28, 'centerDistanceX' : 19.71, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 29, 'centerDistanceX' : 9.85, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6739.22,9759.0,10.54], 'SII6731b' : [6752.67,5619.0,6.687], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 30, 'centerDistanceX' : 19.71, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 31, 'centerDistanceX' : 0.0, 'centerDistanceY' : -5.63, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6690.73,4620.0,9.003], 'SII6731a' : [6704.03,1788.0,3.407], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 32, 'centerDistanceX' : 4.93, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6704.67,5092.0,6.453], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 33, 'centerDistanceX' : 9.85, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6696.25,4962.0,6.091], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 34, 'centerDistanceX' : 14.78, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6691.87,8735.0,8.075], 'SII6731a' : [6706.82,6557.0,7.639], 'SII6716b' : [6733.43,9477.0,8.112], 'SII6731b' : [6748.07,5846.0,6.489], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 35, 'centerDistanceX' : 19.71, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6742.47,5434.0,3.855], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 36, 'centerDistanceX' : -19.71, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.05,2980.0,4.793], 'SII6731a' : [6716.18,1423.0,5.806], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 37, 'centerDistanceX' : -59.12, 'centerDistanceY' : 61.92, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 38, 'centerDistanceX' : 0.0, 'centerDistanceY' : 33.77, 'Halpha6563a': [6549.76,6705.0,3.55], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6736.32,6570.0,4.34], 'SII6731b' : [6750.37,8503.0,6.568], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 39, 'centerDistanceX' : -9.85, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6703.14,2462.0,3.883], 'SII6731a' : [6717.75,1527.0,2.519], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6751.14,999.0,2.391], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 40, 'centerDistanceX' : -4.93, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6701.27,3646.0,3.082], 'SII6731a' : [6716.52,3083.0,5.172], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 41, 'centerDistanceX' : 4.93, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6738.91,2665.0,4.411], 'SII6731b' : [6752.28,3340.0,6.145], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 42, 'centerDistanceX' : 0.0, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6734.28,3962.0,4.035], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 43, 'centerDistanceX' : -4.93, 'centerDistanceY' : 2.81, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 44, 'centerDistanceX' : 0.0, 'centerDistanceY' : 11.26, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [5019.16,22632.0,4.518], 'SII6716a': [6696.69,11508.0,9.678], 'SII6731a' : [6711.16,7869.0,7.848], 'SII6716b' : [6734.63,5445.0,5.453], 'SII6731b' : [6750.64,5318.0,7.499], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 45, 'centerDistanceX' : -9.85, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6699.73,1951.0,3.38], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6742.34,1603.0,4.608], 'SII6731b' : [6752.28,1811.0,3.416], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 46, 'centerDistanceX' : 4.93, 'centerDistanceY' : -2.81, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6694.48,3801.0,8.345], 'SII6731a' : [6708.72,2039.0,5.233], 'SII6716b' : [6726.54,668.0,3.088], 'SII6731b' : [6743.97,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 47, 'centerDistanceX' : 0.0, 'centerDistanceY' : 5.63, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6697.74,3178.0,4.418], 'SII6731a' : [6711.08,1150.0,3.527], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 48, 'centerDistanceX' : 0.0, 'centerDistanceY' : -11.26, 'Halpha6563a': [6549.98,4595.0,2.931], 'Halpha6563b': [6573.68,3883.0,4.169], 'OIII5007a': [4988.7,14171.0,3.14], 'OIII5007b': [5021.75,22884.0,2.656], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6735.77,22561.0,6.976], 'SII6731b' : [6750.14,16430.0,6.323], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 49, 'centerDistanceX' : -4.93, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6692.54,2388.0,3.33], 'SII6731a' : [6705.75,2411.0,3.647], 'SII6716b' : [6730.47,9512.0,7.648], 'SII6731b' : [6745.08,9124.0,8.033], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 50, 'centerDistanceX' : -9.85, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6732.63,1868.0,5.519], 'SII6731b' : [6750.32,3188.0,10.12], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 51, 'centerDistanceX' : -4.93, 'centerDistanceY' : -8.44, 'Halpha6563a': [6558.23,3686.0,2.419], 'Halpha6563b': [6580.97,1491.0,1.793], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6695.66,12368.0,7.381], 'SII6731a' : [6711.25,9037.0,9.048], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6752.33,2033.0,2.843], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 52, 'centerDistanceX' : 0.0, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 53, 'centerDistanceX' : 0.0, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6695.99,7757.0,8.995], 'SII6731a' : [6709.99,5581.0,7.868], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 54, 'centerDistanceX' : -64.05, 'centerDistanceY' : 36.59, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [5003.38,15577.0,2.736], 'OIII5007b': [5010.25,17841.0,3.992], 'SII6716a': [6705.61,5673.0,8.19], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6740.65,5508.0,9.61], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 55, 'centerDistanceX' : -24.64, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6703.13,5475.0,6.71], 'SII6731a' : [6717.19,1449.0,2.726], 'SII6716b' : [6731.1,8508.0,6.565], 'SII6731b' : [6746.4,7373.0,5.942], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 56, 'centerDistanceX' : -9.85, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [6752.45,4574.0,5.247], 'ArIII7136' : [6752.46,3457.0,4.476]},
{'fiber': 57, 'centerDistanceX' : -29.56, 'centerDistanceY' : 33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6699.43,8976.0,12.1], 'SII6731a' : [6717.03,1664.0,2.895], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 58, 'centerDistanceX' : -14.78, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6713.91,4975.0,6.623], 'SII6731a' : [6727.96,5155.0,7.41], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 59, 'centerDistanceX' : -34.49, 'centerDistanceY' : 25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6700.87,6884.0,6.102], 'SII6731a' : [6715.78,3243.0,3.756], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 60, 'centerDistanceX' : -19.71, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6740.37,3235.0,7.354], 'SII6731b' : [6755.33,3839.0,6.6], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 61, 'centerDistanceX' : -24.64, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6698.45,6995.0,5.355], 'SII6731a' : [6711.11,2971.0,3.138], 'SII6716b' : [6736.4,10819.0,6.237], 'SII6731b' : [6750.19,5231.0,4.913], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 62, 'centerDistanceX' : -4.93, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6691.09,6081.0,7.42], 'SII6731a' : [6706.9,3606.0,4.967], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 63, 'centerDistanceX' : -14.78, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6695.05,6801.0,8.332], 'SII6731a' : [6708.12,3949.0,8.927], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 64, 'centerDistanceX' : -14.78, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6696.05,3298.0,7.259], 'SII6731a' : [6714.8,8493.0,12.81], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 65, 'centerDistanceX' : -24.64, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6696.09,4114.0,7.733], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6734.97,4773.0,4.971], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 66, 'centerDistanceX' : -4.93, 'centerDistanceY' : -2.81, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6737.43,7481.0,6.41], 'SII6731b' : [6751.1,963.0,2.34], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 67, 'centerDistanceX' : -14.78, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 68, 'centerDistanceX' : -19.71, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6741.61,3560.0,6.167], 'SII6731b' : [6756.12,3897.0,7.147], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 69, 'centerDistanceX' : -19.71, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 70, 'centerDistanceX' : -64.05, 'centerDistanceY' : -30.96, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6704.1,0.0,0.0], 'SII6731a' : [6719.17,3164.0,5.669], 'SII6716b' : [6728.12,2157.0,5.502], 'SII6731b' : [6741.92,1702.0,3.702], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 71, 'centerDistanceX' : -9.85, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6697.64,6347.0,6.761], 'SII6731a' : [6711.26,2755.0,4.769], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 72, 'centerDistanceX' : -29.56, 'centerDistanceY' : -33.77, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6701.61,2018.0,4.635], 'SII6731a' : [6714.72,1485.0,3.686], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 73, 'centerDistanceX' : -29.56, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6738.43,9211.0,3.504], 'SII6731b' : [6752.36,5277.0,4.732], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 74, 'centerDistanceX' : -34.49, 'centerDistanceY' : -8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6701.5,5182.0,6.961], 'SII6731a' : [6717.77,5797.0,7.725], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 75, 'centerDistanceX' : -29.56, 'centerDistanceY' : 16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.62,2895.0,4.994], 'SII6731a' : [6716.83,5538.0,6.046], 'SII6716b' : [6741.37,722.0,2.374], 'SII6731b' : [6754.7,2541.0,5.616], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 76, 'centerDistanceX' : -34.49, 'centerDistanceY' : 8.44, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6698.82,8039.0,6.714], 'SII6731a' : [6713.97,7082.0,5.429], 'SII6716b' : [6740.03,7786.0,6.309], 'SII6731b' : [6755.38,8260.0,7.74], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 77, 'centerDistanceX' : -29.56, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 78, 'centerDistanceX' : -34.49, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.77,1952.0,2.199], 'SII6731a' : [6709.09,1772.0,2.871], 'SII6716b' : [6737.29,6487.0,4.756], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 79, 'centerDistanceX' : 0.0, 'centerDistanceY' : -16.89, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6702.74,1318.0,3.651], 'SII6731a' : [6718.44,2504.0,3.964], 'SII6716b' : [6738.8,1210.0,2.705], 'SII6731b' : [6753.06,823.0,2.725], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 80, 'centerDistanceX' : -64.05, 'centerDistanceY' : 2.81, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [6705.44,7949.0,7.032], 'SII6731a' : [6721.96,16228.0,8.417], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 81, 'centerDistanceX' : -24.64, 'centerDistanceY' : -25.33, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [6732.72,6155.0,7.314], 'SII6731b' : [6747.07,8696.0,10.28], 'ArIII7136' : [0.0,0.0,0.0]},
{'fiber': 82, 'centerDistanceX' : -19.71, 'centerDistanceY' : 0.0, 'Halpha6563a': [0.0,0.0,0.0], 'Halpha6563b': [0.0,0.0,0.0], 'OIII5007a': [0.0,0.0,0.0], 'OIII5007b': [0.0,0.0,0.0], 'SII6716a': [0.0,0.0,0.0], 'SII6731a' : [0.0,0.0,0.0], 'SII6716b' : [0.0,0.0,0.0], 'SII6731b' : [0.0,0.0,0.0], 'ArIII7136' : [0.0,0.0,0.0]},
]

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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SII6716a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SII6716b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SII6731a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SII6731b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_OIII5007a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_OIII5007b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_Halpha6563a.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_Halpha6563b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SIIa.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_SIIb.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_WIYN_and_GTCa.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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

plt.ylim(-110., 110.)
plt.xlim(-110., 110.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_WIYN_and_GTCb.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
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
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_WIYN_and_GTC_a+b.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()




#sum up spectral range 6690 - 6760 and reproduce an image
wavelengthRange = [6690,6760]
inFile = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMedian.fits'#_6684-6760n.fits'#-skyMean.fits'
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

plt.ylim(-80., 80.)
plt.xlim(-80., 80.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('center distance [arc sec]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_reconstructed_image_WIYN.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

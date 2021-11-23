from myUtils import getImageData,getHeader

fitsFileName = '/Users/azuri/daten/uni/HKU/HASH/DECAPS/imdb.fits'

catHeader = getHeader(fitsFileName)
catData = getImageData(fitsFileName)
print('catData = ',catData)
print('catData = ',catData.size)
print('catHeader = ',len(catHeader),': ',catHeader)
print('type(catData) = ',type(catData))
print('dir(catData) = ',dir(catData))
print('catData.columns = ',catData.columns)
print('catData.getfield("ra") = ',catData.getfield('ra'))

import csvFree,csvData

fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/eugene/hash_FitsFiles_200524.csv')
pngImages = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/eugene/hash_pngimages_200524.csv')

with open('/Users/azuri/daten/uni/HKU/interns_projects/eugene/copy_png_images.sh','w') as f:
    for i in range(fitsFiles.size()):
        if fitsFiles.getData('setname',i) == 'DBS_May2008':
            idPNMain = fitsFiles.getData('idPNMain',i)
            found = False
            print('checking idPNMain ',idPNMain)
            for j in range(pngImages.size()):
                if pngImages.getData('idPNMain',j) == idPNMain:
                    if ('SHS' in pngImages.getData('OUT_DIR',j)) or ('SSS' in pngImages.getData('OUT_DIR',j)):
                        f.write('cp -p '+pngImages.getData('OUT_DIR',j)+pngImages.getData('OutImage',j)+'_overlay.png .\n')
                        found = True
            if not found:
                print('ERROR: could not find an image for idPNMain ',idPNMain)

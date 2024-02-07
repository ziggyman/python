import os
import shutil
import csvFree,csvData

def copyContinuumSources(date):
    areasFNames = [os.path.join(os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/',date),'areas.csv'),
                   os.path.join(os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/',date),'areas.csv')]
    print('areasFNames = ',areasFNames)

    for areasFName in areasFNames:
        print('areasFName = ',areasFName)
        if 'B_data' in areasFName:
            outputDir = os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/CentralStars/B_data/',date)
        else:
            outputDir = os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/CentralStars/R_data/',date)
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        if not os.path.exists(outputDir.replace('B_data','R_data')):
            os.makedirs(outputDir.replace('B_data','R_data'))
        outFileName = os.path.join(outputDir,'to_extract.csv')
        outFileNameExists = os.path.exists(outFileName)
        existingOutFiles = []
        if outFileNameExists:
            existingOutFilesCSV = csvFree.readCSVFile(outFileName)
            existingOutFiles = existingOutFilesCSV.getData('fName')
        with open(outFileName,'a') as f:
            if not outFileNameExists:
                f.write('fName,notes\n')
            areas = csvFree.readCSVFile(areasFName)
            for i in range(areas.size()):
                imageName = areas.getData('fName',i)
                imageName = imageName[:imageName.rfind('.')]+'-sky.fits'
                note = areas.getData('notes',i)
                if imageName[0] != '#':
                    if 'ontinuum' in note:
                        if 'B_data' in imageName:
                            newImageName = os.path.join(outputDir,imageName[imageName.rfind('/')+1:])
                            if newImageName not in existingOutFiles:
                                shutil.copyfile(imageName,newImageName)
                                f.write(newImageName+','+areas.getData('notes',i)+'\n')
                            imageName = imageName.replace('B_data','R_data').replace('b_','r_')
                            newImageName = os.path.join(outputDir.replace('B_data','R_data'),imageName[imageName.rfind('/')+1:])

                            areasB = csvFree.readCSVFile(areasFName.replace('B_data','R_data'))
                            outFileNameB = os.path.join(outputDir.replace('B_data','R_data'),'to_extract.csv')
                            outFileNameBExists = os.path.exists(outFileNameB)
                            existingOutFilesB = []
                            if outFileNameBExists:
                                existingOutFilesBCSV = csvFree.readCSVFile(outFileNameB)
                                existingOutFilesB = existingOutFilesBCSV.getData('fName')
                            for j in range(areasB.size()):
                                if imageName[:imageName.find('dbs')+3] in areasB.getData('fName',j):
                                    imageName = areasB.getData('fName',j)
                                    imageName = imageName[:imageName.rfind('.')]+'-sky.fits'
                                    newImageName = os.path.join(outputDir.replace('B_data','R_data'),imageName[imageName.rfind('/')+1:])
                            with open(outFileNameB,'a') as fB:
                                if not outFileNameBExists:
                                    fB.write('fName,notes\n')
                                if newImageName not in existingOutFilesB:
                                    shutil.copyfile(imageName,newImageName)
                                    fB.write(newImageName+','+areasB.getData('notes',i)+'\n')
                        else:
                            newImageName = os.path.join(outputDir,imageName[imageName.rfind('/')+1:])
                            if newImageName not in existingOutFiles:
                                shutil.copyfile(imageName,newImageName)
                                f.write(newImageName+','+areas.getData('notes',i)+'\n')
                            imageName = imageName.replace('R_data','B_data').replace('r_','b_')
                            newImageName = os.path.join(outputDir.replace('R_data','B_data'),imageName[imageName.rfind('/')+1:])

                            areasB = csvFree.readCSVFile(areasFName.replace('R_data','B_data'))
                            outFileNameB = os.path.join(outputDir.replace('R_data','B_data'),'to_extract.csv')
                            outFileNameBExists = os.path.exists(outFileNameB)
                            existingOutFilesB = []
                            if outFileNameBExists:
                                existingOutFilesBCSV = csvFree.readCSVFile(outFileNameB)
                                existingOutFilesB = existingOutFilesBCSV.getData('fName')
                            for j in range(areasB.size()):
                                if imageName[:imageName.find('dbs')+3] in areasB.getData('fName',j):
                                    imageName = areasB.getData('fName',j)
                                    imageName = imageName[:imageName.rfind('.')]+'-sky.fits'
                                    newImageName = os.path.join(outputDir.replace('R_data','B_data'),imageName[imageName.rfind('/')+1:])
                            with open(outFileNameB,'a') as fB:
                                if not outFileNameBExists:
                                    fB.write('fName,notes\n')
                                if not newImageName in existingOutFilesB:
                                    shutil.copyfile(imageName,newImageName)
                                    fB.write(newImageName+','+areasB.getData('notes',i)+'\n')
if False:
    outputDirB = outFileName[:outFileName.rfind('/')]
    outputDirR = outFileName[:outFileName.rfind('/')].replace('B_data','R_data')
    if not os.path.exists(outputDirB):
        os.makedirs(outputDirB)
    if not os.path.exists(outputDirR):
        os.makedirs(outputDirR)
    outFileNameExists = os.path.exists(outFileName)
    outFileNameRExists = os.path.exists(outFileName.replace('B_data','R_data'))

#        os.remove(outFileName)
    with open(outFileName,'a') as f:
        with open(outFileName.replace('B_data','R_data'),'a') as fR:
            if not outFileNameExists:
                f.write('fName,notes\n')
            if not outFileNameRExists:
                fR.write('fName,notes\n')
            inputList = csvFree.readCSVFile(areasFileNameIn)
            for i in range(inputList.size()):
                imageName = inputList.getData('fName',i)
                imageName = imageName[:imageName.rfind('.')]+'-sky.fits'
                if imageName[0] != '#':
                    note = inputList.getData('notes',i)
                    if 'ontinuum' in note:
                        if inputDir != outputDirB:
                            newImageName = os.path.join(outputDirB,imageName[imageName.rfind('/')+1:])
                            newImageNameR = os.path.join(outputDirR,imageName[imageName.rfind('/')+1:].replace('b_','r_'))
                            shutil.copyfile(imageName,newImageName)
                            f.write(newImageName+','+note+'\n')
                            fR.write(newImageNameR+','+note+'\n')
                            if 'B_data' in imageName:
                                imageName.replace('B_data','R_data')
                                shutil.copyfile(imageName,newImageNameR)
                                imageName = newImageName

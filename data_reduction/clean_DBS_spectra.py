import os

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def findObjectNameInObjectsList(objectName,objectsList):
    for i in range(len(objectsList)):
        if objectName in objectsList[i]:
            print('found ',objectName,' at position ',i)
            return i

path = '/Users/azuri/spectra/MSSSO_2m3_DBS_aug07'
hashList_fileName = os.path.join(path,'hash_red.list')
hashList = readFileToArr(hashList_fileName)
print('hashList = ',hashList)

objectsList_fileName = os.path.join(path,'objects_red.list')
objectsList = readFileToArr(objectsList_fileName)

path = hashList_fileName[:hashList_fileName.rfind('/')]
for line in hashList:
    objectName = line[:line.find('_')]
    print('objectName = ',objectName)
    pos = findObjectNameInObjectsList(objectName,objectsList)

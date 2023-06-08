import csvFree,csvData

def get_IDPNMain_from_name(names,fName_hash_tbCNames,fName_hash_PNMain):
    hash_tbCNames = csvFree.readCSVFile(fName_hash_tbCNames)
    hash_PNMain = csvFree.readCSVFile(fName_hash_PNMain)
    for i in range(hash_tbCNames.size()):
        hash_tbCNames.setData('Name',i,hash_tbCNames.getData('Name',i).lower().replace('_','').replace('j','').replace(' ','').replace('png','').replace('pn','').replace('sp','').replace('mpa','').replace('hen','').replace('he',''))#.replace('-',''))
    newNames = hash_tbCNames.getData('Name')
    with open(fName_hash_tbCNames+'.tmp','w') as f:
        for name in newNames:
            f.write(name+'\n')
#    print('newNames = ',newNames)
    resultIDs = []
    for name in names:
        idPNMain = ''
        idx = hash_tbCNames.find('Name',name.lower().replace('_','').replace('j','').replace(' ','').replace('png','').replace('pn','').replace('sp','').replace('mpa','').replace('hen','').replace('he',''))#.replace('-',''))
#        print('name = ',name,': idx = ',idx)
        if idx[0] < 0:
            print('getIDPNMain_from_name: ERROR: could not find ',name)
#            if name in newNames:
#                print('found <'+name+'> in newNames')
            if name.find('.') != name.rfind('.'):
                if '+' not in name:
                    if '-' not in name:
                        name = name[0:5]+'-'+name[5:]
            idx = hash_PNMain.find('PNG',name)
            if idx[0] < 0:
                print('did not find ',name,' in PNGs either')
                STOP
            idPNMain = hash_PNMain.getData('idPNMain',idx[0])
        else:
            idPNMain = hash_tbCNames.getData('idPNMain',idx[0])
        if len(idx) > 1:
            for i in range(len(idx)-1):
                if hash_tbCNames.getData('idPNMain',idx[i+1]) != idPNMain:
                    print('ERROR: found more than 1 idPNMain: name = ',name,': idPNMain = [',idPNMain,', ',hash_tbCNames.getData('idPNMain',idx[i+1]),']')
                    STOP
        if hash_tbCNames.getData('idPNMain',idx[0]) not in resultIDs:
           resultIDs.append(idPNMain)
        print('name = ',name,': idx = ',idx,', idPNMain = ',idPNMain)
    return resultIDs

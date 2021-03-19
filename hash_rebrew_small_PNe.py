import csvFree
import csvData

diamFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbAngDiam.csv'
pnMainFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain.csv'
outFileName = '/Users/azuri/daten/uni/HKU/HASH/rebrew.sh'

pnMain = csvFree.readCSVFile(pnMainFileName)
diams = csvFree.readCSVFile(diamFileName)

n = 0
with open(outFileName,'w') as f:
    f.write('hashpn brew all ')
    for i in range(pnMain.size()):
        if (pnMain.getData('domain',i) == 'Galaxy') and (pnMain.getData('PNstat',i) in ['T','L','P']):
            idPNMain = pnMain.getData('idPNMain',i)
            for j in range(diams.size()):
                if (diams.getData('idPNMain',j) == idPNMain) and (diams.getData('InUse',j) == '1') and (float(diams.getData('MajDiam',j)) <= 10.):
                    f.write(idPNMain+',')
                    print('i = ',i,': idPNMain = ',idPNMain)
                    n+=1
print('need to rebrew ',n,' PNe')

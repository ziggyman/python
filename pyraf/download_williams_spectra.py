fName = '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/williams/PN_Gallery_Spectra_links.txt'

with open(fName,'r') as f:
    lines = f.readlines()

with open(fName.replace('txt','sh'),'w') as f:
    for line in lines:
        line = line.strip().split()
        f.write('wget '+line[1]+'\n')

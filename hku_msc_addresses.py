import os

#import csvFree,csvData

path = '/Users/azuri/daten/uni/HKU/addresses/'
csvFiles = ['BT Scie Eng.csv','Guangdong SciEng.csv','Nanjing SciEng.csv','ShanghaiSciEng.csv','ShandongSciEng.csv', 'Wuhan SciEng.csv']
for csvFile in csvFiles:
    with open(os.path.join(path,csvFile),'r') as fil:
        lines = fil.readlines()
        outFile = os.path.join(path,csvFile[0:-3]+'txt')

        with open(outFile,'w') as f:
            for line in lines[1:]:
                row = line.split(',')
                if row[0] != '':
                    for txt in row:
                        print('txt = ',txt)
                        f.write(txt+'\n')
                        if txt == '\r':
                            f.write('\n\n\n')

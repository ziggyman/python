def fixTable(fNameIn, fNameOut=None):
    newLines = []
    with open(fNameIn, 'r') as f:
        lines = f.readlines()
        print('lines = ',lines)
        for line in lines:
            line = line.replace(' ± ','±')
            line = line.replace(' ±','±')
            lineA = ''
            lineB = ''
            if line[0:len(line)-5].find('NTT') > 0:
                str = 'NTT'
            else:
                str = 'HST'
            pos = line.find(str)
            lineA = line[0:pos+3]+'\n'
            lineB = line[pos+4:]
            newLines.append(lineA.replace(' ',','))
            newLines.append(lineB.replace(' ',','))
    if fNameOut:
        with open(fNameOut,'w') as f:
            f.write('PNG,Morphology,PA,GPA,Telescope')
            for line in newLines:
                f.write(line)

def crossMatch(fNameInA, fNameInB, fNameOutMatch=None, fNameOutNoMatch=None):
    with open(fNameInA,'r') as fA:
        linesA = fA.readlines()
    with open(fNameInB,'r') as fB:
        linesB = fB.readlines()
    match = []
    noMatch = []
    for lineA in linesA:
        lineATemp = lineA[lineA.find(',')+1:]
        pngNameA = lineATemp[:lineATemp.find(',')]
        print('pngNameA = <'+pngNameA+'>')
        found = False
        newLine = lineA
        for lineB in linesB:
            pngNameB = lineB[:lineB.find(',')]
            print('pngNameB = <'+pngNameB+'>')
            if pngNameA == pngNameB:
                found = True
                newLine = newLine[:newLine.find('\n')] + ',' + lineB
        if found:
            match.append(newLine)
        else:
            noMatch.append(newLine)
    if fNameOutMatch:
        with open(fNameOutMatch,'w') as f:
#            f.write(linesA[0][:linesA[0].find('\n')]+','+linesB[0])
            for line in match:
                f.write(line)
    if fNameOutNoMatch:
        with open(fNameOutNoMatch,'w') as f:
            f.write(linesA[0])
            for line in noMatch:
                f.write(line)

fixTable('/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.txt','/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.csv')
crossMatch('/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe.csv','/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe+Rees_Zijlstra.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe-Rees_Zijlstra.csv')
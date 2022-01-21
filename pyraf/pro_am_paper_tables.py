import csvFree,csvData

csvFiles = ['/Users/azuri/daten/uni/HKU/pro-am-paper/Table_I_210Obj_19122021.csv',
             '/Users/azuri/daten/uni/HKU/pro-am-paper/Table_II_610Obj_19122021.csv',
             '/Users/azuri/daten/uni/HKU/pro-am-paper/Table_III_104Obj_19122021.csv',
             '/Users/azuri/daten/uni/HKU/pro-am-paper/Table_X_Spectra_TLP_221Obj_19122021.csv',
             '/Users/azuri/daten/uni/HKU/pro-am-paper/Table_X+1_Spectra_Mimics_62Obj_19122021.csv',
            ]

for csvFile in csvFiles:
    texFile = csvFile[:-3]+'tex'
    with open(texFile,'w') as f:
        csvData = csvFree.readCSVFile(csvFile)

        oldHeader = csvData.header
        for headerItem in oldHeader:
#            headerItem = headerItem.replace('\r','')
            print('headerItem = ',headerItem)
#            if 'Catalogue' in headerItem:
#                print('headerItem = <'+headerItem+'>')
#                STOP
            if headerItem in ['Mir signal (WISE)',
                                'Name (HASH)',
#                                'Status (HASH)',
                                'MajDiam() (HASH)',
                                'MinDiam() (HASH)',
                                'Catalogue (HASH)\r',
                                'MinDiam (HASH)',
                                'Other name',
                                'Source images'
#                                '',
                                ]:
                csvData.removeColumn(headerItem)
        f.write('\\longtab{\n')
        f.write('\\clearpage\n')
        f.write('\\onecolumn\n')
        f.write('\\begin{landscape}\n')
        if csvFile == csvFiles[0]:
            f.write('\\setlength{\\tabcolsep}{5pt}\n')
        elif csvFile == csvFiles[1]:
            f.write('\\setlength{\\tabcolsep}{5pt}\n')
        elif csvFile == csvFiles[2]:
            f.write('\\setlength{\\tabcolsep}{4pt}\n')
        elif csvFile == csvFiles[3]:
            f.write('\\setlength{\\tabcolsep}{5pt}\n')
        else:
            f.write('\\setlength{\\tabcolsep}{5pt}\n')
        f.write('	\\begin{longtable}{ | *{'+str(len(csvData.header))+')}{l|} }\n')
        if csvFile == csvFiles[0]:
            f.write('	\\caption{Table containing 123 true (T), 51 likely (L) and 36 possible (P) PNe}\n')
            f.write('	\\label{tab1}\\\\\n')
        elif csvFile == csvFiles[1]:
            f.write('	\\caption{Table containing 610 objects awaiting classification}\n')
            f.write('	\\label{tab2}\\\\\n')
        elif csvFile == csvFiles[2]:
            f.write('	\\caption{Table containing 104 objects classified as mimics}\n')
            f.write('	\\label{tab3}\\\\\n')
        elif csvFile == csvFiles[3]:
            f.write('	\\caption{Table of all our combined spectral observations of PN candidates until December 2021}\n')
            f.write('	\\label{tabX}\\\\\n')
        elif csvFile == csvFiles[4]:
            f.write('	\\caption{Table of 66 PN mimics revealed by both professional and amateur spectroscopy of some of our PN candidates}\n')
            f.write('	\\label{tabXI}\\\\\n')
        #f.write('	\\toprule\n')
        f.write('		\\hline\n')
        header = csvData.header
        print('header = ',header)
        for headerItem in header:
            if headerItem == 'Halpha -  [OIII] images  -   French spectrum':
                f.write('\\vtop{\\hbox{\\strut Halpha - [OIII] images -}\\hbox{\\strut French spectrum}}')
            elif headerItem == 'TD* (m)':
                f.write('TelDiam (m)')
            else:
                f.write(headerItem.replace('Dimension','Size').replace('()','(\\arcsec)'))
            if headerItem != header[len(header)-1]:
                f.write(' & ')
        f.write('\\\\\n')
        f.write('		\\endhead  % header material\n')
        f.write('		\\hline\\endfoot  % footer material\n')
        f.write('	\\hline\n')

        for iRow in range(csvData.size()):
            row = csvData.getData(iRow)
            print('row = ',row)
            for iHeaderItem in range(len(header)):
                print('header[iHeaderItem] = ',header[iHeaderItem])
                printRow = row[iHeaderItem].replace('     ',' ').replace('    ',' ').replace('   ',' ').replace('  ',' ').replace('PrivateList','')
                if header[iHeaderItem] == 'Publication':
                    printRow = printRow.replace(' ','; ')
                print('printRow = ',printRow)
                maxColWidth = 20
                #if header[iHeaderItem] == 'Publication':
                #    maxColWidth = 12
                separator = '; '
                if header[iHeaderItem] == 'Halpha -  [OIII] images  -   French spectrum':
                    separator = ' '

                if ((header[iHeaderItem] == 'Source images')
                     or (header[iHeaderItem] == 'Other name') 
                     or (header[iHeaderItem] == 'Publication') 
                     or (header[iHeaderItem] == 'Spectrum') 
                     or (header[iHeaderItem] == 'Halpha -  [OIII] images  -   French spectrum')
                    ) and (len(printRow) > maxColWidth):
                    multirow = []
                    tempRow = printRow.split(separator)
                    print('tempRow = ',tempRow)
                    newRow = tempRow[0]
                    for rowItem in tempRow[1:]:
                        newRow += separator
                        if len(newRow + rowItem) > maxColWidth:
                            multirow.append(newRow)
                            newRow = rowItem
                        else:
                            newRow += rowItem
                        print('newRow = ',newRow)
                    multirow.append(newRow)
                    print('multirow = ',multirow)
                    printRow = '\\vtop{'
                    for cell in multirow:
                        printRow += '\\hbox{\\strut '+cell+'}'
                    printRow += '}'
                    print('printRow = ',printRow)
                    #STOP
                f.write(printRow.replace('_','\_').replace('P. Le D','P. Le D\^u').replace('&','\&').replace('\r','').replace('; -; ','; ').replace(';;',';'))
                if iHeaderItem != len(header)-1:
                    f.write(' & ')
            #STOP
            f.write('\\\\\n')
#            if row[0] == 'SuFe 1':
#                STOP
        f.write('\\\\\n')
        f.write('	\\hline\n')
        f.write('\\end{longtable}\n')
#        if csvFile == csvFiles[0]:
#        f.write('}\n')
        f.write('\\end{landscape}\n')
        f.write('}\n')
#    Stop

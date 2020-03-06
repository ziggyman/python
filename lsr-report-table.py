import numpy as np

googleAnalyticsTableOrig = '/Users/azuri/daten/uni/HKU/website/reports/04Mar2020/users.dat'
latexFile = '/Users/azuri/daten/uni/HKU/website/reports/04Mar2020/users.tex'

with open(googleAnalyticsTableOrig,'r') as f:
    lines = f.readlines()
for i in range(len(lines)):
    lines[i] = lines[i].rstrip('\n')
print(lines)

headerLanguage = lines[0]
print(headerLanguage)

headerCategories = lines[1].split('\t')
print(headerCategories)

header = lines[2].split('\t')
headerColumns = [headerLanguage]
headerColumns.append(header[0])
#headerColumns.append(header[1])
headerColumns.append(header[2])
headerColumns.append(header[4])
print('headerColumns = ',headerColumns)

total = [' ']
total.append(lines[4])
print('total[1] = ',total[1])

#total.append(lines[6])
#print('total[2] = ',total[2])

total.append(lines[8])
print('total[2] = ',total[2])

total.append(lines[12])
print('total[3] = ',total[3])

data = []
print('len(lines) = ',len(lines))
userLanguages = 0
for iLine in np.arange(22,len(lines)-2,10):
    userLanguages += 1
    print('userLanguages = ',userLanguages,': iLine = ',iLine)
    dataRow = [lines[iLine].replace('\t',' ')]
    dataRow.append(lines[iLine+1])
#    dataRow.append(lines[iLine+2])
    dataRow.append(lines[iLine+3])
    dataRow.append(lines[iLine+5])
    data.append(dataRow)
print('data = ',data)

with open(latexFile,'w') as f:
    f.write('\\documentclass{article}[12pt]\n')
    f.write('\\usepackage{graphicx}\n')
    f.write('\\usepackage{tabularx}\n')
    f.write('\\usepackage{vmargin}\n')
    f.write('\\usepackage{hyperref}\n')
    f.write('\\setpapersize[portrait]{A4}\n')
#    f.write('\\setmargins{10mm}{20mm}{16cm}{26cm}{0cm}{0mm}{0mm}{0mm}\n')
    f.write('\\pagestyle{empty}\n')
    f.write('\\begin{document}\n')
    f.write('\\centerline{\Large{Report about the usage of the HKU-LSR website}}\n')
    f.write('\\ \\\\\n')
    f.write('During the day the new \href{https://www.lsr.hku.hk/}{HKU-LSR Website} was activated (May 14th 2019) and today (March 4th 2020), a total of 1,195 different visitors from 44 countries have accessed the website 1,729 times. Thereby each user visited on average 3.1 pages. The number of visitors per day is shown in Fig.~\\ref{fig1} and the details per country in Tab.~\\ref{tab1}.\\\\\n')
    f.write('\\ \\\\\n')
    f.write('\\begin{center}\n')
    f.write('\\begin{figure}[h]\n')
    f.write('\\centering\n')
    f.write('\\includegraphics[width=\\textwidth]{user-numbers}\n')
    f.write('\\caption{Number of visitors per day.\n')
    f.write('}\n')
    f.write('\\label{fig1}\n')
    f.write('\\end{figure}\n')
    f.write('\\begin{table}\n')
    f.write('\\begin{tabularx}{10.2cm}{|c|c|c|c|}\n')
    f.write('\\hline\n')
    f.write(headerColumns[0])
    for col in headerColumns[1:]:
        f.write('&'+col)
    f.write('\\\\\n')
    f.write('(total)')
    for tot in total[1:len(total)-1]:
        f.write('&('+tot+')')
    f.write('& (avg. '+total[len(total)-1]+')')
    f.write('\\\\\n')
    f.write('\\hline\n')
    for lang in data:
        print('lang = ',lang)
        f.write(lang[0])
        for dat in lang[1:]:
            dat = dat.replace('%','\\%')
            print('dat = ',dat)
            f.write('&'+dat)
        f.write('\\\\\n')
    f.write('\\hline\n')

#    f.write('\\n')

    f.write('\\end{tabularx}\n')
    f.write('\\caption{Usage details per country of origin of the visitors.\n')
    f.write('}\n')
    f.write('\\label{tab1}\n')
    f.write('\\end{table}\n')
    f.write('\\end{center}\n')
    f.write('\\end{document}\n')

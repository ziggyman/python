import csv

tableName = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH-True-PN-EPA-measures.csv'

table = csv.DictReader(open(tableName))

ones = []
twos = []
threes = []
for row in table:
    difference = row['ZIG-QAP(EPA)']
    if (difference == '#VALUE!') or (difference == '') or (difference == ' '):
        difference = '0'
    #print('difference = ',difference)
    difference = float(difference)
    flag = int(row['flag'])
    id = row['HASH ID']
    if (abs(difference) > 50) and (abs(180.-abs(difference)) > 50.):
        print('id = ',id,': EPA Ziggy = ',row['EPA Ziggy'],', EPA QAP = ',row['QAP EPA'],': difference = ',difference,': quality flag = ',flag)
        if flag == 1:
            ones.append(id)
        elif flag == 2:
            twos.append(id)
        else:
            threes.append(id)

print('found ',len(ones),' with flag 1')
print('found ',len(twos),' with flag 2')
print('found ',len(threes),' with flag 3')

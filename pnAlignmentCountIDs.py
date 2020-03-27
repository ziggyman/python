import os

path = '/Users/azuri/daten/uni/HKU/PN alignment'
outf = os.path.join(path,'uniqueIDs_B2.csv')
with open(outf,'r') as tf:
    line = tf.readlines()[0].split(',')
b2 = [int(i) for i in line]
print('b2 = ',b2)

outf = os.path.join(path,'uniqueIDs_E2.csv')
with open(outf,'r') as tf:
    line = tf.readlines()[0].split(',')
e2 = [int(i) for i in line]
print('e2 = ',e2)

outf = os.path.join(path,'uniqueIDs_BE2.csv')
with open(outf,'r') as tf:
    line = tf.readlines()[0].split(',')
be2 = [int(i) for i in line]
print('be2 = ',be2)

for i in b2:
    if i not in be2:
        print(i,' found in b2 but not found in be2')
    if i in e2:
        print(i,' found in b2 and in e2')
for i in e2:
    if i not in be2:
        print(i,' found in e2 but not found in be2')
    if i in b2:
        print(i,' found in e2 and in b2')

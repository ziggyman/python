import datetime
import csvFree,csvData

fileName = '/Users/azuri/daten/uni/HKU/HASH/analytics/hash_userlist.csv'

userdata = csvFree.readCSVFile(fileName)
nActiveLastYear = 0
# date in yyyy/mm/dd format
today = datetime.datetime(2022, 8, 29)
oneYearAgo = datetime.datetime(2021, 8, 29)
for i in range(userdata.size()):
    lastLog = userdata.getData('lastLog',i)
    if (lastLog != 'NULL') and ('0000' not in lastLog):
        year,month,day = lastLog.split(' ')[0].split('-')
        print('i = ',i,': year = ',year,', month = ',month,', day = ',day)
        thisDate = datetime.datetime(int(year), int(month), int(day))
        if thisDate > oneYearAgo:
            nActiveLastYear += 1
print('found ',nActiveLastYear,' users active within the last 12 months')

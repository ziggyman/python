from datetime import datetime,time
from matplotlib import pyplot as plt

import csvFree
import csvData

csvFileIn = '/Users/azuri/daten/uni/HKU/website/reports/04Mar2020/hash-userlog.csv'
userList = '/Users/azuri/daten/uni/HKU/website/reports/04Mar2020/hash-userlist.csv'
startDate = '2019-03-04'

data = csvFree.readCSVFile(csvFileIn,',',False)
print('header = ',data.header)
print('data.size() = ',data.size())

userCSV = csvFree.readCSVFile(userList,',',False)

start = False
users = []
dates = []
times = []
for i in range(data.size()):
    if (not start) and (data.getData('time',i).split(' ')[0] == startDate):
        start = True
    if start:
        user = data.getData('userName',i)
        if user not in users:
            users.append(user)

        dateAndTime = data.getData('time',i)
        date,timeStr = dateAndTime.split(' ')
        year, month, day = date.split('-')
        dates.append(datetime(int(year), int(month), int(day)))
        hour, min, sec = timeStr.split(":")
        print('hour = ',hour)
        print('min = ',min)
        print('sec = ',sec)
        times.append(time(int(hour), int(min), int(sec)))
print('len(users) = ',len(users))
#print('dates = ',dates[0:10])
print('times = ',times[0:10])

affiliations = []
for user in users:
    for i in range(userCSV.size()):
        if user == userCSV.getData('userName',i):
            affiliation = userCSV.getData('affiliation',i)
            if affiliation not in affiliations:
                affiliations.append(affiliation)
print('len(affiliations) = ',len(affiliations))

plt.figure(figsize=(12, 6))
plt.hist(dates, bins=365)#, rwidth=0.8)
plt.xlabel('Date')
plt.ylabel('Number of queries per day')
plt.show()

#plt.hist(times, bins=1440)#, rwidth=0.8)
#plt.xlabel('Time')
#plt.ylabel('Number of queries pe')
#plt.gcf().autofmt_xdate()
#plt.show()
hour_list = [t.hour for t in times]
print(hour_list)
numbers=[x for x in range(0,24)]
labels=map(lambda x: str(x), numbers)
plt.xticks(numbers, labels)
plt.xlim(0,23)
plt.hist(hour_list,bins=24)#,align='left')
plt.xlabel('Hong Kong time')
plt.ylabel('Number of queries for the last 12 months')
plt.xticks(range(24),['0:00','1:00','2:00','3:00','4:00','5:00','6:00','7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'])
plt.show()


import datetime
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
import numpy as np
import csvFree,csvData

fileName = '/Users/azuri/daten/uni/HKU/HASH/analytics/2023/hash_userlist.csv'

userdata = csvFree.readCSVFile(fileName)
nActiveLastYear = 0
# date in yyyy/mm/dd format
today = datetime.datetime(2023, 8, 29)
oneYearAgo = datetime.datetime(2022, 9, 1)
for i in range(userdata.size()):
    lastLog = userdata.getData('lastLog',i)
    if (lastLog != 'NULL') and ('0000' not in lastLog):
        year,month,day = lastLog.split(' ')[0].split('-')
        print('i = ',i,': year = ',year,', month = ',month,', day = ',day)
        thisDate = datetime.datetime(int(year), int(month), int(day))
        if thisDate > oneYearAgo:
            nActiveLastYear += 1
print('found ',nActiveLastYear,' users active within the last 12 months out of ',userdata.size(),' registered users')

g3s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/data-export_G3-3years.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/usersG3-4years.csv')]
g4s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/data-export_G4-2.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/usersG4-2.csv')]

countries_g3s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/countries-G3.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/countries-G3.csv')]
countries_g4s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/countries-G4.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/countries-G4.csv')]

startDates = [datetime.date(2020,9,1),datetime.date(2019,9,1)]
startDatesG4 = [datetime.date(2023,4,1),datetime.date(2023,7,1)]
for k in range(2):
    g3 = g3s[k]
    g4 = g4s[k]
    countries_g3 = countries_g3s[k]
    countries_g4 = countries_g4s[k]

    delta = datetime.timedelta(days=1)
    startDate = startDates[k]
    midDate = datetime.date(2022,8,31)
    endDate = datetime.date(2023,8,31)
    delta = endDate - startDate
    datesOld = []
    for i in range(delta.days+1):
        datesOld.append(startDate + datetime.timedelta(days=i))
    print('datesOld = ',datesOld)
    delta = endDate - midDate
    datesNew = []
    for i in range(delta.days+1):
        datesNew.append(midDate + datetime.timedelta(days=i))
    print('datesNew = ',datesNew)
    nUsersOld = np.zeros(len(datesOld))
    nUsersNew = np.zeros(len(datesNew))
    for i in range(g3.size()):
        date = g3.getData('Day Index',i)
        month = int(date[:date.find('/')])
        day = int(date[date.find('/')+1:date.rfind('/')])
        year = 2000+int(date[date.rfind('/')+1:])
        date = datetime.date(year,month,day)
        print('date = ',date,': day = ',day,', month = ',month,', year = ',year)
        found = False
        for j in range(len(datesOld)):
            if datesOld[j] == date:
                nUsersOld[j] = int(g3.getData('1 Day Active Users',i))
                found = True
        for j in range(len(datesNew)):
            if datesNew[j] == date:
                nUsersNew[j] = int(g3.getData('1 Day Active Users',i))
                found = True
        if not found:
            print('ERROR: date ',date,' not found')


    print('g4.header = ',g4.header)
    startDateG4 = startDatesG4[k]
    endDate = datetime.date(2023,8,31)
    delta = endDate - startDateG4
    datesG4 = []
    for i in range(delta.days+1):
        datesG4.append(startDateG4 + datetime.timedelta(days=i))
    print('datesG4 = ',datesG4)
    nUsersG4 = np.zeros(len(datesG4))
    for i in range(g4.size()):
        date = startDateG4 + datetime.timedelta(days=int(g4.getData('Nth day',i)))
        for j in range(len(datesOld)):
            if datesOld[j] == date:
                nUsersOld[j] = int(g4.getData('Users\r',i).strip('\r'))
                found = True
        for j in range(len(datesNew)):
            if date == datesNew[j]:
                nUsersNew[j] = int(g4.getData('Users\r',i).strip('\r'))

    #idx = np.argsort(dates)
    print('len(datesNew) = ',len(datesNew),', len(datesOld) = ',len(datesOld))

    allUsers = np.sum(np.array(nUsersNew))
    print('all together ',allUsers,' users visited last year')
    #plt.plot(datesNew,nUsersNew,'b-',label='01/09/2022 to 31/08/2023')
    plt.plot(datesOld,nUsersOld)#,'g-',label='01/09/2021 to 31/08/2022')
    #plt.plot(datesG4,nUsersG4)#,'c-',label='G4')
    #plt.legend()
    plt.xlim(datesOld[0],datesOld[len(datesOld)-1])
    plt.ylabel('Users per day')
    #datesP = matplotlib.dates.date2num(dates)
    #plt.plot_date(datesP, nUsers)
    plt.show()

    countries = []
    users = []
    sessions = []
    for j in range(25):
        countries.append(countries_g3.getData('Country',j))
        users.append(int(countries_g3.getData('Users',j)))
        sessions.append(int(countries_g3.getData('Sessions',j)))
        for l in range(countries_g4.size()):
            country = countries_g4.getData('Country',l)
            if country == countries[len(countries)-1]:
                users[len(users)-1] += int(countries_g4.getData('Users',l))
                sessions[len(sessions)-1] += int(countries_g4.getData('Engaged sessions',l))
    print('users = ',users)
    idx = np.array([-user for user in users]).argsort()
    print("idx = ",idx)
    countries = np.array(countries)[idx]
    users = np.array(users)[idx]
    sessions = np.array(sessions)[idx]
    print('Country, Users, Sessions')
    for j in range(21):
        print(countries[j],', ',users[j],', ',sessions[j])

import datetime
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
import pandas as pd
import numpy as np
import csvFree,csvData

fileName = '/Users/azuri/daten/uni/HKU/HASH/analytics/2024/hash_userlist.csv'

userdata = csvFree.readCSVFile(fileName)
nActiveLastYear = 0
# date in yyyy/mm/dd format
today = datetime.datetime(2023, 8, 29)
oneYearAgo = datetime.datetime(2022, 9, 1)
for i in range(userdata.size()):
    lastLog = userdata.getData('lastLog',i)
    if (lastLog != 'NULL') and ('0000' not in lastLog):
        year,month,day = lastLog.split(' ')[0].split('-')
        #print('i = ',i,': year = ',year,', month = ',month,', day = ',day)
        thisDate = datetime.datetime(int(year), int(month), int(day))
        if thisDate > oneYearAgo:
            nActiveLastYear += 1
print('found ',nActiveLastYear,' users active within the last 12 months out of ',userdata.size(),' registered users')

g3s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/data-export_G3-3years.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/usersG3-4years.csv')]
g4s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/data-export_G4-2.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/usersG4-2.csv')]

users2024 = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2024/users.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2024/users.csv')]

countries_g3s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/countries-G3.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/countries-G3.csv')]
countries_g4s = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2023/countries-G4.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2023/countries-G4.csv')]
countries2024 = [csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/analytics/2024/countries.csv'),csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/website/analytics/2024/countries.csv')]


startDates = [datetime.date(2020,9,1),datetime.date(2019,9,1)]
startDatesG4 = [datetime.date(2023,4,1),datetime.date(2023,7,1)]
startDatesLatest = [datetime.date(2023,9,1),datetime.date(2023,9,1)]
for k in range(2):
    g3 = g3s[k]
    g4 = g4s[k]
    users = users2024[k]
    countriesLatest = countries2024[k]

    delta = datetime.timedelta(days=1)
    startDate = startDates[k]
    midDate = datetime.date(2022,8,31)
    endDate = datetime.date(2023,8,31)
    newEndDate = datetime.date(2024,8,31)

    delta = endDate - startDate
    datesOld = []
    for i in range(delta.days+1):
        datesOld.append(startDate + datetime.timedelta(days=i))
    print('k = ',k,': datesOld = ',datesOld)

    delta = endDate - midDate
    datesNew = []
    for i in range(delta.days+1):
        datesNew.append(midDate + datetime.timedelta(days=i))
#    print('k = ',k,': datesNew = ',datesNew)

    delta = newEndDate - startDate
    allDates = []
    for i in range(delta.days+1):
        allDates.append(startDate + datetime.timedelta(days=i))

    startDateLatest = startDatesLatest[k]
    endDateLatest = datetime.date(2024,8,31)
    delta = endDateLatest - startDateLatest
    datesLatest = []
    for i in range(delta.days+1):
        datesLatest.append(startDateLatest + datetime.timedelta(days=i))
#    print('k = ',k,': datesLatest = ',datesLatest)

    nUsersOld = np.zeros(len(datesOld))
    nUsersNew = np.zeros(len(datesNew))
    nUsersLatest = np.zeros(len(datesLatest))
    allUsers = np.zeros(len(allDates))
    for i in range(g3.size()):
        date = g3.getData('Day Index',i)
        month = int(date[:date.find('/')])
        day = int(date[date.find('/')+1:date.rfind('/')])
        year = 2000+int(date[date.rfind('/')+1:])
        date = datetime.date(year,month,day)
        print('k = ',k,': date = ',date,': day = ',day,', month = ',month,', year = ',year)
        found = False
        for j in range(len(datesOld)):
            if datesOld[j] == date:
                nUsersOld[j] = int(g3.getData('1 Day Active Users',i))
                found = True
        for j in range(len(datesNew)):
            if datesNew[j] == date:
                nUsersNew[j] = int(g3.getData('1 Day Active Users',i))
                found = True
        for j in range(len(allDates)):
            if allDates[j] == date:
                allUsers[j] = int(g3.getData('1 Day Active Users',i))
        if not found:
            print('k = ',k,': ERROR: date ',date,' not found')


    startDateG4 = startDatesG4[k]
    endDate = datetime.date(2023,8,31)
    delta = endDate - startDateG4
    datesG4 = []
    for i in range(delta.days+1):
        datesG4.append(startDateG4 + datetime.timedelta(days=i))
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
        for j in range(len(datesLatest)):
            if date == datesLatest[j]:
                nUsersLatest[j] = int(g4.getData('Users\r',i).strip('\r'))
        for j in range(len(allDates)):
            if date == allDates[j]:
                allUsers[j] = int(g4.getData('Users\r',i).strip('\r'))

    delta = endDateLatest - startDateLatest
    for i in range(users.size()):
        date = startDateLatest + datetime.timedelta(days=int(users.getData('Nth day',i)))
        for j in range(len(datesLatest)):
            if datesLatest[j] == date:
                nUsersLatest[j] = int(users.getData('Users\r',i))
                found = True
        for j in range(len(allDates)):
            if allDates[j] == date:
                allUsers[j] = int(users.getData('Users\r',i))

    print('k = ',k,': all together ',np.sum(np.array(nUsersLatest)),' users visited last year')
#    print('datesNew = ',datesNew)
    datesNew = [dateNew + datetime.timedelta(days=365) for dateNew in datesNew]
    mydpi = 150
    plt.figure(figsize=(1500/mydpi,600/mydpi),dpi=mydpi)
    plt.plot(datesNew,pd.DataFrame(nUsersNew).rolling(min_periods=1,window=7,center=True).mean().values,'b-',label='01/09/2022 to 31/08/2023')
    plt.plot(datesLatest,pd.DataFrame(nUsersLatest).rolling(min_periods=1,window=7,center=True).mean().values,'g-',label='01/09/2023 to 31/08/2024')
    #plt.plot(datesG4,nUsersG4)#,'c-',label='G4')
    plt.legend()
    plt.xlim(datesLatest[0],datesLatest[len(datesLatest)-1])
    plt.ylabel('Users per day')
    #datesP = matplotlib.dates.date2num(dates)
    #plt.plot_date(datesP, nUsers)
    plt.tight_layout()
    if k == 0:
        plt.savefig('/Users/azuri/daten/uni/HKU/HASH/analytics/2024/users_per_day_HASH_last_2_years_comparison.png',dpi=mydpi)
    else:
        plt.savefig('/Users/azuri/daten/uni/HKU/website/analytics/2024/users_per_day_LSR_last_2_years_comparison.png',dpi=mydpi)
    plt.show()

    mydpi = 150
    plt.figure(figsize=(1500/mydpi,600/mydpi),dpi=mydpi)
    rollingAv = pd.DataFrame(allUsers).rolling(min_periods=1, window=7, center=True).mean()
    print('rollingAv = ',type(rollingAv),': ',rollingAv)
    print('dir(rollingAv) = ',dir(rollingAv))
    plt.plot(allDates,rollingAv.values,'b-',label='01/09/2022 to 31/08/2023')
    #plt.plot(datesLatest,nUsersLatest,'g-',label='01/09/2023 to 31/08/2024')
    #plt.plot(datesG4,nUsersG4)#,'c-',label='G4')
    #plt.legend()
    plt.xlim(startDate,allDates[len(allDates)-1])
    plt.ylabel('Users per day')
    #datesP = matplotlib.dates.date2num(dates)
    #plt.plot_date(datesP, nUsers)
    plt.tight_layout()
    if k == 0:
        plt.savefig('/Users/azuri/daten/uni/HKU/HASH/analytics/2024/users_per_day_HASH_full.png',dpi=mydpi)
    else:
        plt.savefig('/Users/azuri/daten/uni/HKU/website/analytics/2024/users_per_day_LSR_full.png',dpi=mydpi)
    plt.show()

    countries = []
    users = []
    sessions = []
    for j in range(countriesLatest.size()):
        countries.append(countriesLatest.getData('Country',j))
        users.append(int(countriesLatest.getData('Users',j)))
        sessions.append(int(countriesLatest.getData('Engaged sessions',j)))
    print('k = ',k,': users = ',np.sum(users))
    if k == 0:
        idx = np.array([-user for user in users]).argsort()
    else:
        idx = np.array([-session for session in sessions]).argsort()
    print("k = ',k,': idx = ",idx)
    countries = np.array(countries)[idx]
    users = np.array(users)[idx]
    sessions = np.array(sessions)[idx]
    print('k = ',k,': Country, Users, Sessions')
    for j in range(21):
        if k == 0:
            print(countries[j],', ',users[j],', ',sessions[j])
        else:
            print(countries[j],', ',sessions[j],', ',users[j])

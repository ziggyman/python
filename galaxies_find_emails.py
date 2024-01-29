import csvFree,csvData

iaus323 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/email_list_2nd announcement_IAUS323.txt','\t',False)
for i in range(iaus323.size()):
    iaus323.setData('email\r',i,iaus323.getData('email\r',i).replace('\r','').replace(' ',''))

apnVII_speakers = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/APNVII-invited-speaker-emails.txt')
#first_names_apnVII_speakers = apnVII_speakers.getData('first name')
#last_names_apnVII_speakers = apnVII_speakers.getData('last name')
#emails_apnVII_speakers = apnVII_speakers.getData('email')

apnVII_participants = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/APNVII-reg-participants-240817.list')
#first_names_apnVII = apnVII_participants.getData('first name')
#names_apnVII = apnVII_participants.getData('last name')

apnV = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/APNV-participants.txt')

apnVI = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/APNVI_participants.csv','\t',False)

workplansII = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/emails_WorkPLANsII.csv')

q_pn_distribution_list = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/PN-email-distribution-list copy.csv')
for i in range(q_pn_distribution_list.size()):
    for key in q_pn_distribution_list.header:
        q_pn_distribution_list.setData(key,i,q_pn_distribution_list.getData(key,i).strip())
csvFree.writeCSVFile(q_pn_distribution_list,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/PN-email-distribution-list copy.csv')

all_previous = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/cfp-mailing-list.csv')
print('all_previous = ',all_previous.header)
#names_all_previous = all_previous.getData('Name')
#first_names_all_previous = [name[name.find(';')+1:].strip() for name in names_all_previous]
#names_all_previous = [name[:name.find(';')] for name in names_all_previous]
#all_previous.removeColumn('last name')
#all_previous.removeColumn('first name')
#all_previous.addColumn('last name',names_all_previous)
#all_previous.addColumn('first name',first_names_all_previous)
#csvFree.writeCSVFile(all_previous,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/cfp-mailing-list.csv')
#emails_all_previous = all_previous.getData('EMail')
#print('names_all_previous = ',names_all_previous)
#print('emails_all_previous = ',emails_all_previous)

all = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/CFP-GE_Oct09.txt','\t',False)
#print('all = ',all.header)
#print(all.data)
#names_all = all.getData('Name')
#first_names_all = [name[name.find(',')+1:].strip() for name in names_all]
#names_all = [name[:name.find(',')] for name in names_all]
#emails_all = all.getData('EMail')
#print('names_all = ',names_all)
#print('emails_all = ',emails_all)
#all.removeColumn('first name')
#all.removeColumn('last name')
#all.addColumn('first name',first_names_all)
#all.addColumn('last name',names_all)
#csvFree.writeCSVFile(all,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/CFP-GE_Oct09.txt','\t')

cfp1023 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/cfp-GE 10.23.csv','\t',False)
#print('cfp1023 = ',cfp1023.header)
#print(cfp1023.data)
#names_cfp1023 = cfp1023.getData('Name')
#first_names_cfp1023 = [name[name.find(';')+1:].strip() for name in names_cfp1023]
#names_cfp1023 = [name[:name.find(';')] for name in names_cfp1023]
#emails_cfp1023 = cfp1023.getData('EMail')
#print('names_cfp1023 = ',names_cfp1023)
#print('emails_cfp1023 = ',emails_cfp1023)
#cfp1023.removeColumn('first name')
#cfp1023.removeColumn('last name')
#cfp1023.addColumn('first name',first_names_cfp1023)
#cfp1023.addColumn('last name',names_cfp1023)
#csvFree.writeCSVFile(cfp1023,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/cfp-GE 10.23.csv','\t')

new = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/possible_authors.csv')
#print('new = ',new.header)
#print(new.data)
#names_new = new.getData('last name')
#first_names_new = new.getData('first name')
#print('names_new = ',names_new)
#emails_new = new.getData('email\r')
#print('emails_new = ',emails_new)
#emails_new = [email.replace('\r','') for email in emails_new]

iaus384 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/emails_IAUS384.csv')
#names_iaus384 = iaus384.getData('Surname')
#first_names_iaus384 = iaus384.getData('Name')
#emails_iaus384 = iaus384.getData('Email\r')
#print('names_iaus384 = ',names_iaus384)
#print('emails_iaus384 = ',emails_iaus384)

iaus283 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/iaus283_participants.txt','\t',False)
#names_iaus283 = iaus283.getData('name')
#first_names_iaus283 = [name[name.find(',')+2:] for name in names_iaus283]
#names_iaus283 = [name[0:name.find(',')] for name in names_iaus283]
#print('names_iaus283 = <'+names_iaus283[0]+'>, ',names_iaus283)

apn8 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/apn8.csv')
#print('apn8 = ',apn8.header)
#names_apn8 = apn8.getData('last name')
#first_names_apn8 = apn8.getData('first name')
#print('names_apn8 = ',names_apn8)

nature = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/nature_paper_authors.csv')
#names_nature = nature.getData('last name')
#first_names_nature = nature.getData('first name')
#print('names_nature = ',names_nature)


hash = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/hash_users.csv')

#names_iaus323 = iaus323.getData('Last Name')
#first_names_iaus323 = iaus323.getData('First Name')
#emails_iaus323 = iaus323.getData('email\r')
#print('names_iaus323 = ',names_iaus323)
#print('fist_names_iaus323 = ',first_names_iaus323)
#print('emails_iaus323 = ',emails_iaus323)

final_list = csvData.CSVData()
final_list.header = ['first name','last name','email']
for i in range(iaus384.size()):
    final_list.append([iaus384.getData('Name',i),iaus384.getData('Surname',i),iaus384.getData('Email\r',i).strip()])
print('with IAUS384: final_list.size() = ',final_list.size())

for i in range(apnVII_speakers.size()):
    if final_list.find('email',apnVII_speakers.getData('email',i))[0] < 0:
        final_list.append([apnVII_speakers.getData('first name',i),apnVII_speakers.getData('last name',i),apnVII_speakers.getData('email',i).strip()])
print('with apnVII_speakers: final_list.size() = ',final_list.size())

for i in range(apnVI.size()):
    if final_list.find('email',apnVI.getData('email',i))[0] < 0:
        final_list.append([apnVI.getData('First name',i),apnVI.getData('Last name',i),apnVI.getData('email',i).strip()])
print('with apnVI: final_list.size() = ',final_list.size())

for i in range(new.size()):
    if final_list.find('email',new.getData('email\r',i))[0] < 0:
        final_list.append([new.getData('first name',i),new.getData('last name',i),new.getData('email\r',i).strip()])
print('with possible authors: final_list.size() = ',final_list.size())

for i in range(workplansII.size()):
    if final_list.find('email',workplansII.getData('email\r',i))[0] < 0:
        final_list.append([workplansII.getData('first name',i),workplansII.getData('last name',i),workplansII.getData('email\r',i).strip()])
print('with WorkPLANsII: final_list.size() = ',final_list.size())

for i in range(q_pn_distribution_list.size()):
    if final_list.find('email',q_pn_distribution_list.getData('email',i))[0] < 0:
        final_list.append([q_pn_distribution_list.getData('first name',i),q_pn_distribution_list.getData('last name',i),q_pn_distribution_list.getData('email',i).strip()])
print('with q_pn_distribution_list: final_list.size() = ',final_list.size())

for i in range(iaus323.size()):
    if final_list.find('email',iaus323.getData('email\r',i))[0] < 0:
        final_list.append([iaus323.getData('First Name',i),iaus323.getData('Last Name',i),iaus323.getData('email\r',i).strip()])
print('with IAUS323: final_list.size() = ',final_list.size())

for i in range(hash.size()):
    if final_list.find('email',hash.getData('email',i))[0] < 0:
        if hash.getData('firstName',i) == 'NULL':
            final_list.append(['HASH user','',hash.getData('email',i)])
        else:
            name = hash.getData('firstName',i)
            final_list.append([name[:name.find(' ')],name[name.find(' ')+1:],hash.getData('email',i)])

print('with HASH users: final_list.size() = ',final_list.size())

csvFree.writeCSVFile(final_list,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/final_list.csv')

not_found = []

#apnVII_speakers add names
if False:
    names_and_emails_lists = [[q_pn_distribution_list,'email','first name','last name'],
                            [all,'EMail','fist name','last name'],
                            [all_previous,'EMail','fist name','last name'],
                            [cfp1023,'EMail','first name','last name']]
    apnVII_speakers.addColumn('first name',['' for i in range(apnVII_speakers.size())])
    apnVII_speakers.addColumn('last name',['' for i in range(apnVII_speakers.size())])
    for find_names_list in names_and_emails_lists:
        csv,emailKey,firstNameKey,lastNameKey = find_names_list
        for i in range(apnVII_speakers.size()):
            email = apnVII_speakers.getData('email',i)
            found = csv.find(emailKey,email)[0]
            if found >= 0:
                print('found email <'+email+'> in csv at position ',found)
                if len(apnVII_speakers.getData('first name',i)) < len(csv.getData(firstNameKey,found)):
                    apnVII_speakers.setData('first name',i,csv.getData(firstNameKey,found))
                if len(apnVII_speakers.getData('last name',i)) < len(csv.getData(lastNameKey,found)):
                    apnVII_speakers.setData('last name',i,csv.getData(lastNameKey,found))
        for i in range(csv.size()):
            if final_list.find('email',csv.getData(emailKey,i))[0] < 0:
                final_list.append([csv.getData(firstNameKey,i),csv.getData(lastNameKey,i),csv.getData(emailKey,i)])
    csvFree.writeCSVFile(apnVII_speakers,'/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/APNVII-invited-speaker-emails.txt')

    STOP

names_to_find = [names_apn8,names_iaus283,names_apnVII,names_nature]
first_names_to_find = [first_names_apn8, first_names_iaus283, first_names_nature]
for j in range(len(names_to_find)):
    name_to_find = names_to_find[j]
    first_name_to_find = first_names_to_find[j]
    for k in range(len(name_to_find)):
        name_apn8 = name_to_find[k]
        first_name = first_name_to_find[k]
        print('searching for name <'+name_apn8+'>')
        found = False
        for i in range(len(names_all)):
#            if 'Sabin' in names_all[i]:
#                print('found Sabin in names_all[',i,'] = ',names_all[i])
            if (name_apn8 in names_all[i]):# or (names_all[i] in name_apn8):
                print('found name <'+name_apn8+'> in names_all[',i,'] = ',names_all[i],', email address = <'+emails_all[i]+'>')
                found = True
                if name_apn8 not in final_list.getData('last name'):
                    final_list.append([first_names_all[i],name_apn8,emails_all[i]])
                else:
                    print('name <'+name_apn8+'> found in names_all but already in final list')
        for i in range(len(names_all_previous)):
#            if 'Sabin' in names_all_previous[i]:
#                print('found Sabin in names_all_previous[',i,'] = ',names_all_previous[i])
            if ((name_apn8 in names_all_previous[i])# or (names_all_previous[i] in name_apn8)
                ) and not found:
                print('found name <'+name_apn8+'> in names_all_previous[',i,'] = ',names_all_previous[i],', email address = <'+emails_all_previous[i]+'>')
                found = True
                if name_apn8 not in final_list.getData('last name'):
                    final_list.append([first_names_all_previous[i],name_apn8,emails_all_previous[i]])
                else:
                    print('name <'+name_apn8+'> found in names_all_previous but already in final list')
        for i in range(len(names_iaus384)):
#            if 'Sabin' in names_iaus384[i]:
#                print('found Sabin in names_iaus384[',i,'] = ',names_iaus384[i])
            if ((name_apn8 in names_iaus384[i])# or (names_iaus384[i] in name_apn8)
                ) and not found:
                print('found name <'+name_apn8+'> in names_iaus384[',i,'] = ',names_iaus384[i],', email address = <'+emails_iaus384[i]+'>')
                found = True
                if name_apn8 not in final_list.getData('last name'):
                    final_list.append([first_names_iaus384[i],name_apn8,emails_iaus384[i]])
                else:
                    print('name <'+name_apn8+'> found in names_iaus384 but already in final list')
        for i in range(len(names_new)):
#            if 'Sabin' in names_new[i]:
#                print('found Sabin in names_new[',i,'] = ',names_new[i])
            if ((name_apn8 in names_iaus384[i])# or (names_iaus384[i] in name_apn8)
                ) and not found:
                print('found name <'+name_apn8+'> in names_new[',i,'] = ',names_new[i],', email address = <'+emails_new[i]+'>')
                found = True
                if name_apn8 not in final_list.getData('last name'):
                    final_list.append([first_names_new[i],name_apn8,emails_new[i]])
                else:
                    print('name <'+name_apn8+'> found in names_new but already in in final list')
        if not found:
            print('did not find email address for name <'+name_apn8+'>')
            if name_apn8 not in not_found:
                not_found.append(first_name.strip(' ')+' '+name_apn8)
            if name_apn8 == 'Sabin':
                STOP

print('final_list.size() = ',final_list.size())
print('not_found: ',len(not_found),': ',not_found)

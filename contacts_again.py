import numpy as np
import csvFree,csvData

eventRegistrationFile = '/Users/azuri/daten/parties/Wild Woods/Event Registration Wild Woods.csv'
oldContactsFile = '/Users/azuri/daten/parties/Wild Woods/contacts/contacts.csv'
contactsToImport = '/Users/azuri/daten/parties/Wild Woods/contacts/contacts_to_import.csv'
existingContactsFile = '/Users/azuri/daten/parties/Wild Woods/existing_contacts.csv'

csvEventRegistration = csvFree.readCSVFile(eventRegistrationFile)
csvOldContacts = csvFree.readCSVFile(oldContactsFile)
existingContacts = csvFree.readCSVFile(existingContactsFile)

csvContactsToImport = csvData.CSVData()
csvContactsToImport.header = csvOldContacts.header

newLine = ['' for i in csvOldContacts.header]

for i in range(csvEventRegistration.size()):
    csvContactsToImport.append(newLine)
    csvContactsToImport.setData('Name',i,csvEventRegistration.getData('Name',i)+' Wild Woods Guest')
    csvContactsToImport.setData('Group Membership',i,'* myContacts')
    csvContactsToImport.setData('Phone 1 - Type',i,'Mobile')
    csvContactsToImport.setData('Phone 1 - Value',i,csvEventRegistration.getData('WhatsApp number - for payment payment instructions; directions will be sent out via WhatsApp broadcast so make sure you add the number +852 5963 3756 to your contacts',i))

for i in np.arange(csvContactsToImport.size()-1,-1,-1):
    print('checking contact ',i)
    if existingContacts.find('Phone 1 - Value',csvContactsToImport.getData('Phone 1 - Value',i),0)[0] > -1:
        print('found ',csvContactsToImport.getData('Name',i))
        csvContactsToImport.removeRow(i)


csvFree.writeCSVFile(csvContactsToImport,contactsToImport)

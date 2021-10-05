import numpy as np
import csvFree,csvData

eventRegistrationFile = '/Users/azuri/daten/parties/Wild Wood/contacts/Completed orders export.csv'
oldContactsFile = '/Users/azuri/daten/parties/Wild Wood/contacts/contacts.csv'
contactsToImport = '/Users/azuri/daten/parties/Wild Wood/contacts/contacts_to_import.csv'
existingContactsFile = '/Users/azuri/daten/parties/Wild Wood/existing_contacts.csv'
volunteersFile = '/Users/azuri/daten/parties/Wild Wood/Wild Woods Volunteers Questionnaire.csv'

csvEventRegistration = csvFree.readCSVFile(eventRegistrationFile)
csvOldContacts = csvFree.readCSVFile(oldContactsFile)
existingContacts = csvFree.readCSVFile(existingContactsFile)
volunteers = csvFree.readCSVFile(volunteersFile)

csvContactsToImport = csvData.CSVData()
csvContactsToImport.header = csvOldContacts.header

newLine = ['' for i in csvOldContacts.header]

print('csvEventRegistration.header = ',csvEventRegistration.header)

iContact = 0
names = []
for i in range(csvEventRegistration.size()):
    name = csvEventRegistration.getData('First Name (Billing)',i)+' '+csvEventRegistration.getData('Last Name (Billing)',i)+' Wild Wood Guest'
    print('name = ',name)
    if name not in names:
        print(name+' not found')
        names.append(name)
        csvContactsToImport.append(newLine)
        csvContactsToImport.setData('Name',iContact,name+' Wild Wood ' + 'Guest')
        csvContactsToImport.setData('Group Membership',iContact,'* myContacts')
        csvContactsToImport.setData('Phone 1 - Type',iContact,'Mobile')
        csvContactsToImport.setData('Phone 1 - Value',iContact,csvEventRegistration.getData('Phone (Billing)',i))
        iContact += 1
        print('csvContactsToImport.size() = ',csvContactsToImport.size())
    else:
        print(name+' found')
        print('csvContactsToImport.size() = ',csvContactsToImport.size())

#for i in np.arange(csvContactsToImport.size()-1,-1,-1):
#    print('checking contact ',i)
#    if existingContacts.find('Phone 1 - Value',csvContactsToImport.getData('Phone 1 - Value',i),0)[0] > -1:
#        print('found ',csvContactsToImport.getData('Name',i))
#        csvContactsToImport.removeRow(i)


csvFree.writeCSVFile(csvContactsToImport,contactsToImport)

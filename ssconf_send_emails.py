import smtplib,ssl
#import email
from email.message import EmailMessage
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
import numpy as np
import csvFree,csvData

emailList = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/conferences/SS2Conf/emails - Sheet1_new.csv')
participantsList = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/conferences/SS2Conf/Registration Form _Space Sustainability Conference, 2024(asia_sat).csv')
path_to_pdf = '/Users/azuri/daten/uni/HKU/conferences/SS2Conf/logo.jpg'

emails_participants = participantsList.getData('email')
emails_general = emailList.getData('email')

n_sent = 0

if True:
    for i in range(1):#participantsList.size()):
        if True:#participantsList.getData('email',i).strip() not in emails_general[0:39]:
            firstName = participantsList.getData('Name',i).strip()
            lastName = participantsList.getData('Surname',i).strip()
            to = participantsList.getData('email',i).strip()

            print('firstName = <'+firstName+'>')
            print('to = <'+to+'>')
        #    STOP
        #firstName = 'Ziggy'

            if firstName == 'Madam/':
                firstName = 'potential'
                lastName = 'delegate'

            email_text = f"""
        Dear %s %s,

        Please find attached the conference logo you will need to show to campus security in order to enter HKU.

        We look forward to welcoming you to Hong Kong and our wonderful Hong Kong University Campus.

        Warm regards,

        Quentin (HKU) and Jean-Paul (EPFL)
        SOC co-chairs
            """ % (firstName,lastName)

            MAIL_USERNAME = "ssclsr"
            MAIL_APP_PASSWORD = "pCdDrDe8pc"

            recipients = [to]
        #    ccs = ["rachel.shang@mdpi.com"]
            print('sending mail to firstName = ',firstName,': to = <'+recipients[0]+'>')
            print('email text = <'+email_text+'>')

            msg = MIMEMultipart()
            msg.attach(MIMEText(email_text, "plain"))
            #msg = MIMEText(email_text)
            msg["Subject"] = 'International Conference on Space Sustainability 2024 at HKU'
            msg["To"] = ", ".join(recipients)
            msg["From"] = f"{MAIL_USERNAME}@hku.hk"
            with open(path_to_pdf,'rb') as f:
                attach = MIMEApplication(f.read(),_subtype="pdf")
            attach.add_header('Content-Disposition','attachment',filename=str(path_to_pdf[path_to_pdf.rfind('/')+1:]))
            msg.attach(attach)
            context = ssl.create_default_context()
            smtp_server = smtplib.SMTP('smtproam.hku.hk', 587)#465)
        #    smtp_server = smtplib.SMTP_SSL('smtproam.hku.hk', 587)#465)
            smtp_server.ehlo()
            smtp_server.starttls(context=context)
            smtp_server.ehlo()
            smtp_server.login(MAIL_USERNAME, MAIL_APP_PASSWORD)
            smtp_server.sendmail(msg["From"], msg["To"], msg.as_string())#["From"], recipients + ccs, msg.as_string())
        #    smtp_server.send_message(msg)#["From"], recipients + ccs, msg.as_string())
            smtp_server.quit()
            n_sent += 1
        else:
            print('already sent email to ',to)

if False:
    for i in np.arange(60,emailList.size(),1):
        if emailList.getData('email',i) not in emails_participants:
            firstName = emailList.getData('Name',i).strip()
            lastName = emailList.getData('Surname',i).strip()
            to = emailList.getData('email',i).strip()

            print('firstName = <'+firstName+'>')
            print('to = <'+to+'>')
        #    STOP
        #firstName = 'Ziggy'

            if firstName == 'Madam/':
                firstName = 'potential'
                lastName = 'delegate'

            email_text = f"""
        Dear %s %s,

        Please find attached the final announcement for the upcoming Space Debris and Sustainability conference Dec 2-4th, 2024, in the vibrant harbour city of Hong Kong : https://ssconf.space/

        The conference is now less than a week away and the program has been finalised. It promises to be an interesting and impactful meeting judging by the list of attendees and talks.

        All talk slots are filled, however late registrations can still be accepted.

        We look forward to welcoming you to Hong Kong and our wonderful Hong Kong University Campus.

        Warm regards,

        Quentin (HKU) and Jean-Paul (EPFL)
        SOC co-chairs
            """ % (firstName,lastName)

            MAIL_USERNAME = "ssclsr"
            MAIL_APP_PASSWORD = "pCdDrDe8pc"

            recipients = [to]
        #    ccs = ["rachel.shang@mdpi.com"]
            print('sending mail to firstName = ',firstName,': to = <'+recipients[0]+'>')
            print('email text = <'+email_text+'>')

            msg = MIMEMultipart()
            msg.attach(MIMEText(email_text, "plain"))
            #msg = MIMEText(email_text)
            msg["Subject"] = 'International Conference on Space Sustainability 2024 at HKU'
            msg["To"] = ", ".join(recipients)
            msg["From"] = f"{MAIL_USERNAME}@hku.hk"
            with open(path_to_pdf,'rb') as f:
                attach = MIMEApplication(f.read(),_subtype="pdf")
            attach.add_header('Content-Disposition','attachment',filename=str(path_to_pdf[path_to_pdf.rfind('/')+1:]))
            msg.attach(attach)
            context = ssl.create_default_context()
            smtp_server = smtplib.SMTP('smtproam.hku.hk', 587)#465)
        #    smtp_server = smtplib.SMTP_SSL('smtproam.hku.hk', 587)#465)
            smtp_server.ehlo()
            smtp_server.starttls(context=context)
            smtp_server.ehlo()
            smtp_server.login(MAIL_USERNAME, MAIL_APP_PASSWORD)
            smtp_server.sendmail(msg["From"], msg["To"], msg.as_string())#["From"], recipients + ccs, msg.as_string())
        #    smtp_server.send_message(msg)#["From"], recipients + ccs, msg.as_string())
            smtp_server.quit()
            n_sent += 1
        else:
            print('already emailed to ',emailList.getData('email',i))


print('sent ',n_sent,' emails. emailList.size() = ',emailList.size())

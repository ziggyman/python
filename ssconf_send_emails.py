import smtplib,ssl
#import email
from email.message import EmailMessage
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart

import csvFree,csvData

emailList = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/conferences/SS2Conf/emails - Sheet1_new.csv')
path_to_pdf = '/Users/azuri/daten/uni/HKU/conferences/SS2Conf/Space-Sustainability-2nd-announcement-final.pdf'

n_sent = 0
for i in range(emailList.size()):
    firstName = emailList.getData('Name',i).strip()
    lastName = emailList.getData('Surname',i).strip()
    to = emailList.getData('email',i).strip()

    print('firstName = <'+firstName+'>')
    print('to = <'+to+'>')
#    STOP
#firstName = 'Ziggy'

    email_text = f"""
Dear %s %s,

Please find enclosed our second announcement for our International Conference on “Space Sustainability”. Please note early-bird registrations (that save you $100USD for standard registrations) are set to close at the end of August.

We encourage early registrations  given the meeting is capped at 120 participants.
The conference will take place from December 2nd to 4th, 2024 on the main campus of HKU with EPFL as co-host.

Please see our conference website for updated details: https://ssconf.space/

This conference is of significant interest where Hong Kong’s special status as an East-West super connector should lead to interesting opportunities for discussions and interactions across the broad areas of concern around Space Debris.

Thank you very much for your attention, and we hope to see you at the conference!

This message  is sent on behalf of the co-chairs of the Scientific Organising Committee, Prof. Quentin Parker from the Laboratory for Space Research at the University of Hong Kong (HKU) and Prof. Jean-Paul Kneib from the Laboratory of Astrophysics and Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland


Laboratory for Space Research
Faculty of Science
The University of Hong Kong
Email: ssclsr@hku.hk
Website: https://ssconf.space/
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
    msg["Subject"] = 'Invitation to the International Conference on Space Sustainability 2024 at HKU'
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
print('sent ',n_sent,' emails. emailList.size() = ',emailList.size())
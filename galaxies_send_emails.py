import smtplib
#import email
from email.message import EmailMessage
from email.mime.text import MIMEText

import csvFree,csvData

emailList = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/final_list_with_most_names.csv')

for i in range(emailList.size()):
    firstName = emailList.getData('first name',i)
    to = emailList.getData('email',i)

#firstName = 'Ziggy'

    email_text = f"""
    Dear %s,



    First I would like to wish you a very successful, happy, and healthy new year 2024.


    I am writing to inform you about the Special Issue “Origins and Models of Planetary Nebulae” which will appear in the journal Galaxies.


    Special issue link: https://www.mdpi.com/journal/galaxies/special_issues/2TVP474A4Q


    As a respected scholar in the field, we would like to extend an invitation to you to submit a research paper or comprehensive review paper to this Special Issue.


    The Special Issue aims to showcase the latest advancements and innovative research in Planetary Nebulae. Your contribution will help make this a truly comprehensive and valuable resource for the research community.


    Please note that there is a high chance that publishing fees will be waived for your contributing paper, however the decision on the waiver will only be made after your paper has been submitted.



    Here is a simple guide to paper submission, for your convenience:

    1. First-time users are required to register before making submissions at https://susy.mdpi.com, you can find the instructions at https://www.mdpi.com/journal/galaxies/instructions.

    2. Enter your account, click "Submit Manuscript" under Submissions Menu.

    3. Input manuscript details from Steps 1 to 5:

    Journal: Galaxies

    Research Topics: Origins and Models of Planetary Nebulae

    4. You will receive an auto-reply letter confirming your successful submission.



    Please let me know if you have any questions or if there is anything else I can assist you with. I look forward to the possibility of working with you on this exciting project.



    With warmest regards,

    Dr. Andreas Ritter, Guest Editor, Laboratory for Space Research, Faculty of Science, The University of Hong Kong

    Prof. Xuan Fang, Guest Editor, National Astronomical Observatories, Chinese Academy of Sciences (NAOC)


    """ % (firstName)

    GMAIL_USERNAME = "azuri.ritter"
    GMAIL_APP_PASSWORD = "qxqvtozcafsmasoh"

    recipients = [to]
    ccs = ["rachel.shang@mdpi.com"]
    print('sending mail to firstName = ',firstName,': to = <'+recipients[0]+', cc = ',ccs[0],'>')
    print('email text = <'+email_text+'>')

    msg = EmailMessage()
    msg.set_content(email_text)
    #msg = MIMEText(email_text)
    msg["Subject"] = 'Feature Paper Invitation [Galaxies] (IF 2.5, ISSN 2075-4434) — Special Issue " Origins and Models of Planetary Nebulae "'
    msg["To"] = ", ".join(recipients)
    msg["From"] = f"{GMAIL_USERNAME}@gmail.com"
    msg["Cc"] = ", ".join(ccs)

    smtp_server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    smtp_server.login(GMAIL_USERNAME, GMAIL_APP_PASSWORD)
    smtp_server.send_message(msg)#["From"], recipients + ccs, msg.as_string())
    smtp_server.quit()

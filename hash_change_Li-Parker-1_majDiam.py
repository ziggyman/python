import csvFree,csvData

csvList = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/Li-Parker-1.csv')

x = 1
cron = 16838
with open('/Users/azuri/daten/uni/HKU/HASH/Li-Parker-1.sql','w') as f:
    for i in range(csvList.size()):
        if csvList.getData('Catalogue',i) == 'Li-Parker-1':
            f.write("INSERT INTO `MainPNData`.`LiParker1`(`idLiParker1`,`idPNMain`,`mapFlag`) VALUES (%d,%d,'y');\n" % (x,int(csvList.getData('idPNMain',i))))
            x+=1
            f.write("UPDATE `MainGPN`.`tbAngDiam` SET `MajDiam`=60 WHERE `idPNMain`=%d;\n" % (int(csvList.getData('idPNMain',i))))
            f.write("INSERT INTO `MainPNUsers`.`cronJobs`(`idcronJobs`,`user`,`idPNMain`,`cronScript`,`parameters`,`date_subm`,`priority`,`zombie`) VALUES (%d,'ziggy',%d,'fetch all %d','a:2:{s:1:\"w\";s:5:\"force\";s:3:\"vvv\";s:0:\"\";}','2024-05-16 15:13:04',99,'n');\n" % (cron,int(csvList.getData('idPNMain',i)),int(csvList.getData('idPNMain',i))))
            cron+=1
            f.write("INSERT INTO `MainPNUsers`.`cronJobs`(`idcronJobs`,`user`,`idPNMain`,`cronScript`,`parameters`,`date_subm`,`priority`,`zombie`) VALUES (%d,'ziggy',%d,'brew all %d','a:2:{s:1:\"w\";s:0:\"\";s:3:\"vvv\";s:0:\"\";}','2024-05-16 15:13:04',99,'n');\n" % (cron,int(csvList.getData('idPNMain',i)),int(csvList.getData('idPNMain',i))))
            cron+=1

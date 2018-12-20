execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate

tab = [
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ194645', 'SkyLeft':[1110,1155], 'SkyRight':[1351,1371], 'ObjectAreas':[[1157,1240],[1262,1330]]},#very wide, find sky from other imeage nearby
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0002/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-15,1.3e-15], 'OName': 'Kn19', 'SkyLeft':[1141,1166], 'SkyRight':[1474,1504], 'ObjectAreas':[[1200,1228],[1251,1306],[1321,1340],[1358,1407],[1425,1450]]},#looks good
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0003/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-15,0.2e-15], 'OName': 'Kn45', 'SkyLeft':[1184,1228], 'SkyRight':[1320,1374], 'ObjectAreas':[[1238,1296]]},#okay, diffuse extended object with continuum
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0004/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'DSHJ205811', 'SkyLeft':[889,991], 'SkyRight':[1570,1608], 'ObjectAreas':[[1090,1112],[1122,1256],[1326,1524],[1540,1557]]},#faint line at 5005, very extended object with sky from both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0005/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-14,0.2e-14], 'OName': 'DSHJ221013', 'SkyLeft':[680,780], 'SkyRight':[1626,1704], 'ObjectAreas':[[1090,1292],[1334,1380],[1404,1429],[1479,1506],[1538,1554]]},#okay, line at 6562, very extended object with sky from both CCDs, possible extended cetral source at x=1455?
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0006/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.15e-15], 'OName': 'IPHASXJ194745', 'SkyLeft':[1246,1264], 'SkyRight':[1428,1443], 'ObjectAreas':[[1300,1396]]},#hmmm, run that and try with sky from other image - did try that, still nothing to be seen here
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0007/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-14,0.1e-14], 'OName': 'Ou3', 'SkyLeft':[1123,1132], 'SkyRight':[1503,1515], 'ObjectAreas':[[1133,1150],[1189,1226],[1241,1300],[1318,1331],[1342,1389],[1407,1432],[1445,1468],[1494,1502]]},#looks good, line at 5004, very extended
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0008/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-14,0.4e-14], 'OName': 'Kn24', 'SkyLeft':[885,918], 'SkyRight':[1753,1790], 'ObjectAreas':[[940,1038],[1090,1180],[1198,1303],[1320,1556],[1566,1700]]},#looks good, line at 5004, very extended, on both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0009/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-16,2.0e-16], 'OName': 'LDu1', 'SkyLeft':[984,1025], 'SkyRight':[1837,1918], 'ObjectAreas':[[1106,1134],[1150,1194],[1254,1320],[1338,1348],[1496,1594]]},#looks good, line at 5004 and 6582, very extended, sky from both CCDs, residual sky lines in upper part
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0010/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-15,0.5e-15], 'OName': 'IPHASXJ193617', 'SkyLeft':[1263,1288], 'SkyRight':[1423,1479], 'ObjectAreas':[[1293,1307],[1320,1331]]},#looks good, pretty compact
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0011/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-15,0.5e-15], 'OName': 'IPHASXJ194728', 'SkyLeft':[1177,1197], 'SkyRight':[1395,1414], 'ObjectAreas':[[1248,1273],[1287,1306],[1323,1350],[1365,1374]]},#looks good, lines at 6716,...
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0012/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-14,0.8e-14], 'OName': 'IPHASXJ194301', 'SkyLeft':[1218,1248], 'SkyRight':[1333,1363], 'ObjectAreas':[[1272,1325]]},#looks good, pretty compact, lines at 6585, 6562
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0012a/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.5e-15], 'OName': 'IPHASXJ194301', 'SkyLeft':[1244,1269], 'SkyRight':[1357,1385], 'ObjectAreas':[[1284,1304],[1324,1348]]},#looks good, pretty compact, lines at 6585, 6562
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'PNG124.1-01.9', 'SkyLeft':[1150,1177], 'SkyRight':[1376,1399], 'ObjectAreas':[[1239,1280],[1299,1350]]},#looks good, lines at 4859 and 6716, strong cosmic ray at 3868, actually quite a few residual cosmics, check image
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0002/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-14,0.8e-14], 'OName': 'IPHASXJ014238+600947', 'SkyLeft':[785,816], 'SkyRight':[1930,1956], 'ObjectAreas':[[864,953],[968,1038],[1091,1126],[1138,1303],[1318,1424],[1440,1558],[1564,1595],[1612,1680],[1703,1720],[1759,1902]]},#looks good, lines at 6584, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0003/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ031058.8+624755', 'SkyLeft':[243,324], 'SkyRight':[1961,1985], 'ObjectAreas':[[363,457],[478,542],[559,576],[594,726],[743,786],[834,855],[872,897],[910,1038],[1090,1194],[1221,1273],[1328,1436],[1505,1560],[1596,1760],[1785,1850],[1877,1909],[1928,1950]]},#looks good at 4957 and 6579, very extended object over both CCDs, residual sky lines in the far red
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0004/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-15,1.5e-15], 'OName': 'IPHASXj040721.5+512422', 'SkyLeft':[424,522], 'SkyRight':[1515,1567], 'ObjectAreas':[[1197,1303],[1320,1450]]},#looks good, lines at 4959, 5004, 6580, 6716, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0005/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-13,0.1e-13], 'OName': 'IPHASXJ051152.2+302751', 'SkyLeft':[662,700], 'SkyRight':[1884,1925], 'ObjectAreas':[[701,1035],[1124,1154],[1171,1300],[1325,1555],[1607,1880]]},#faint emission at 6562, 6580, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0005a/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ051152.2+302751', 'SkyLeft':[662,700], 'SkyRight':[1005,1038], 'ObjectAreas':[[701,1000]]},#laft side of CCD, faint emission at 6562, 6580, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0005b/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ051152.2+302751', 'SkyLeft':[1500,1555], 'SkyRight':[1884,1925], 'ObjectAreas':[[1607,1880]]},#right side of CCD, faint emission at 6562, 6580, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0006/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.15e-14,0.1e-14], 'OName': 'IPHASXJ053650.8+245616', 'SkyLeft':[1170,1199], 'SkyRight':[1511,1541], 'ObjectAreas':[[1201,1300],[1320,1474]]},#faint lines around 6716
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0009/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-14,1.5e-14], 'OName': 'IPHASXJ055242.8+262116', 'SkyLeft':[1233,1248], 'SkyRight':[1534,1545], 'ObjectAreas':[[1265,1294],[1327,1466]]},#looks good, compact, 4959, 5008, 6563
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0013/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-2.0e-16,1.5e-16], 'OName': 'IPHASXJ185322.1+083018', 'SkyLeft':[1213,1236], 'SkyRight':[1542,1557], 'ObjectAreas':[[1256,1283],[1322,1358],[1414,1448],[1491,1537]]},#shell visible at 5005
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0020/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-14,0.5e-14], 'OName': 'IPHASXJ191306.1+025248', 'SkyLeft':[728,743], 'SkyRight':[1836,1878], 'ObjectAreas':[[754,801],[846,870],[882,977],[985,1002],[1011,1038],[1132,1183],[1196,1223],[1245,1304],[1346,1390],[1519,1534],[1565,1626],[1635,1696],[1714,1746],[1758,1806]]},#clear shell at 6583, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0026/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-15,1.2e-15], 'OName': 'SNR060707', 'SkyLeft':[768,828], 'SkyRight':[1698,1742], 'ObjectAreas':[[888,958],[1101,1152],[1180,1210],[1245,1299],[1309,1350],[1368,1399],[1425,1478],[1491,1560]]},#lines at 6561, 6584, without extended regions near the center
##{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0026a/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-1.0e-15,1.2e-15], 'OName': 'SNR060707', 'SkyLeft':[325,513], 'SkyRight':[1936,2021], 'ObjectAreas':[[888,958],[1101,1299],[1309,1350],[1368,1399],[1425,1478],[1491,1560]]},#lines at 6561, 6584, with extended regions near the center
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0027/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-15,0.5e-15], 'OName': 'SNR060715', 'SkyLeft':[1149,1179], 'SkyRight':[1720,1750], 'ObjectAreas':[[1256,1300],[1325,1388],[1566,1634],[1665,1701]]},#ok, some lines at 6563, 6717, bright object bleeding into the nebula, mysterious
#{'FName': 'IPHAS_GTC_DATA/DATA3_2017B/GTC8-17BMEX/OB0028/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,0.5e-15], 'OName': 'SNR060930', 'SkyLeft':[1002,1035], 'SkyRight':[1559,1589], 'ObjectAreas':[[1126,1242],[1258,1348],[1408,1555]]},#possible line at 6566, mysterious
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'KK26', 'SkyLeft':[1092,1109], 'SkyRight':[1568,1602], 'ObjectAreas':[[1112,1283],[1305,1559]]},#very nice blobs, actual lines at cut out areas at 6298 and 6362, maybe use spectrum  without cut-out sky lines
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0002/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-15,0.5e-15], 'OName': 'Kn34', 'SkyLeft':[1110,1160], 'SkyRight':[1450,1484], 'ObjectAreas':[[1163,1183],[1218,1234],[1270,1295],[1325,1357],[1394,1443]]},#ok, lines at 4861, 4960, 5006, problematic sky lines at 6300 and 6330
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0003/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-15,1.2e-15], 'OName': 'CR1', 'SkyLeft':[872,919], 'SkyRight':[1926,2002], 'ObjectAreas':[[1109,1158],[1167,1179],[1190,1263],[1284,1320],[1361,1456],[1470,1513],[1522,1529],[1543,1561],[1564,1625]]},#ok, lines at 4959, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0004/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-17,5.0e-17], 'OName': 'Kn11', 'SkyLeft':[1253,1277], 'SkyRight':[1345,1377], 'ObjectAreas':[[1278,1300],[1318,1340]]},#okay, compact, lines 6547, 6562, 6582
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0005/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'HaWe3', 'SkyLeft':[1140,1225], 'SkyRight':[1395,1422], 'ObjectAreas':[[1229,1298],[1323,1391]]},#okay, lines at 4958, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0006/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Kn48', 'SkyLeft':[1120,1190], 'SkyRight':[1442,1470], 'ObjectAreas':[[1216,1253],[1265,1301],[1321,1332],[1348,1389],[1410,1420]]},#okay lines at 4861, 4958, 5006
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0007/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-2.0e-15,9.0e-15], 'OName': 'Kn58', 'SkyLeft':[1150,1190], 'SkyRight':[1434,1461], 'ObjectAreas':[[1195,1300],[1324,1337],[1358,1430]]},#okay, lines at 4955, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0008/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-15,1.9e-15], 'OName': 'Ou2', 'SkyLeft':[1114,1133], 'SkyRight':[1528,1552], 'ObjectAreas':[[1133,1346],[1405,1443],[1460,1475],[1494,1507]]},#okay, lines at 4955, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0009/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Abell68', 'SkyLeft':[1153,1168], 'SkyRight':[1432,1477], 'ObjectAreas':[[1188,1305],[1317,1368],[1379,1397]]},#okay, lines at 4684, 4861, 4958, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0010/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-14,0.5e-14], 'OName': 'Kn24', 'SkyLeft':[860,939], 'SkyRight':[1703,1758], 'ObjectAreas':[[945,998],[1013,1036],[1090,1241],[1254,1303],[1338,1400],[1453,1558],[1581,1702]]},#hmmm, possible line at 6585, very extended object over both CCDs
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0011/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ190333_Green_Slit', 'SkyLeft':[718,748], 'SkyRight':[1825,1842], 'ObjectAreas':[[1164,1212],[1241,1252],[1268,1282],[1307,1353],[1367,1407],[1435,1445],[1495,1516],[1530,1556],[1602,1650],[1670,1700]]},#faint line at 5005, 6565, 6585, 6829, 6862
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0012/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.21e-14,0.4e-14], 'OName': 'K1-6', 'SkyLeft':[490,521], 'SkyRight':[1773,1824], 'ObjectAreas':[[523,741],[756,784],[799,871],[910,927],[944,1036],[1090,1132],[1148,1233],[1398,1582]]},#okay, lines at 4861, 4959, 5005, very extended object over both CCDs, plus bleeding bright star
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0012a/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'peculiar object', 'SkyLeft':[1091,1133], 'SkyRight':[1149,1186], 'ObjectAreas':[[1134,1148]]},#peculiar object
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0013/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.6e-15,1.8e-15], 'OName': 'LDu1', 'SkyLeft':[495,632], 'SkyRight':[1753,1792], 'ObjectAreas':[[1090,1166],[1182,1226],[1238,1304],[1315,1350],[1395,1452],[1465,1480],[1492,1550],[1564,1573],[1588,1642]]},#very extended object, need to take sky from the left CCD, bright object in left CCD bleeding into nebula/sky
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0015/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,0.2e-15], 'OName': 'IPHASXJ195627', 'SkyLeft':[1200,1239], 'SkyRight':[1340,1373], 'ObjectAreas':[[1259,1337]]},#Hmmm, faint line at 6563
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0016/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-16,2.7e-16], 'OName': 'IRAS20084', 'SkyLeft':[1128,1206], 'SkyRight':[1681,1780], 'ObjectAreas':[[1254,1270],[1290,1299],[1325,1329]]},#nothing really but HASH shows compact cloud in 432
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0017/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,0.7e-15], 'OName': 'IPHASXJ201058', 'SkyLeft':[1213,1234], 'SkyRight':[1370,1396], 'ObjectAreas':[[1288,1340]]},#faint lines at 6565 and 6585
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0018/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,2.0e-15], 'OName': 'Kn36', 'SkyLeft':[366,547], 'SkyRight':[1884,1940], 'ObjectAreas':[[1188,1300],[1345,1361],[1398,1420]]},#ok lines at 4956, 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0019/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'We2-260', 'SkyLeft':[681,760], 'SkyRight':[1782,1901], 'ObjectAreas':[[1090,1240],[1262,1716]]},#line at 6582[800,957],[978,1038],
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0020/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ225420', 'SkyLeft':[1145,1180], 'SkyRight':[1440,1470], 'ObjectAreas':[[1192,1232],[1269,1432]]},#line at 6559, 6580
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0023/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Kn25', 'SkyLeft':[1096,1119], 'SkyRight':[1511,1553], 'ObjectAreas':[[1171,1210],[1265,1280],[1334,1343],[1358,1370],[1398,1468]]},#lines at 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0024/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'We1-10', 'SkyLeft':[415,488], 'SkyRight':[1928,1965], 'ObjectAreas':[[752,839],[860,978],[994,1020],[1090,1195],[1209,1295],[1325,1402],[1423,1472],[1504,1554],[1570,1618],[1633,1658],[1679,1845]]},# line at 4959
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0025/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,4.2e-15], 'OName': 'KTC1', 'SkyLeft':[1144,1190], 'SkyRight':[1390,1455], 'ObjectAreas':[[1244,1301],[1326,1375]]},#easy, lines at 4957 and 5005
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0026/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-15,2.0e-15], 'OName': 'Kn30', 'SkyLeft':[1178,1225], 'SkyRight':[1372,1394], 'ObjectAreas':[[1301,1326],[1340,1354]]},#easy, lines at 4859, 4957, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0027/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Kn31', 'SkyLeft':[734,900], 'SkyRight':[1611,1720], 'ObjectAreas':[[1167,1222],[1240,1301],[1322,1495]]},#okay, lines at 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0028/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': '', 'SkyLeft':[688,956], 'SkyRight':[1764,1878], 'ObjectAreas':[[974,1038],[1128,1188],[1206,1288],[1308,1514],[1537,1554],[1570,1662]]},#lines at 5005, 6560, 6580
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0029/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-14,0.08e-14], 'OName': 'IPHASXJ230323', 'SkyLeft':[802,1031], 'SkyRight':[1696,1866], 'ObjectAreas':[[1150,1222],[1251,1276],[1301,1397]]},#faint line at 6561
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0030/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-15,1.5e-15], 'OName': 'IPHASXJ010133', 'SkyLeft':[520,716], 'SkyRight':[1685,1794], 'ObjectAreas':[[1259,1375]]},#line at 4859
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0031/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.8e-15,3.8e-15], 'OName': 'Ju1', 'SkyLeft':[185,390], 'SkyRight':[1819,1983], 'ObjectAreas':[[790,820],[830,849],[861,871],[901,974],[1004,1035],[1101,1146],[1168,1201],[1222,1236],[1256,1301],[1344,1410],[1432,1495],[1515,1543],[1728,1756],[1772,1789]]},#line at 1960
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0032/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.4e-15,1.0e-15], 'OName': 'Kn51', 'SkyLeft':[391,688], 'SkyRight':[1847,1978], 'ObjectAreas':[[1200,1300],[1319,1448]]},#line at 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0033/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-16,1.5e-16], 'OName': 'FsMv1', 'SkyLeft':[695,817], 'SkyRight':[1650,1809], 'ObjectAreas':[[1131,1164],[1176,1186]]},#difficult, possible line at 6581, without extended objects
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0033a/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-2.0e-16,2.0e-16], 'OName': 'FsMv1', 'SkyLeft':[695,817], 'SkyRight':[1650,1809], 'ObjectAreas':[[801,891],[1131,1164]]},#difficult, possible line at 6581, with left extended object
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0033b/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-0.4e-15,0.7e-15], 'OName': 'FsMv1', 'SkyLeft':[695,817], 'SkyRight':[1650,1809], 'ObjectAreas':[[1131,1164],[1205,1306]]},#difficult, possible line at 6581, with right extended object
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0034/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-16,1.5e-16], 'OName': 'KnJ0240', 'SkyLeft':[1142,1244], 'SkyRight':[1395,1530], 'ObjectAreas':[[1284,1298],[1329,1351]]},#compact, lines at 4988, 5090, 5137
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0035/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,4.0e-16], 'OName': 'KnJ1857', 'SkyLeft':[810,913], 'SkyRight':[1874,2014], 'ObjectAreas':[[1221,1288],[1337,1374]]},#can't see lines just a possibly extended (compact) object in the center
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0036/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.8e-15], 'OName': 'IPHASXJ214032', 'SkyLeft':[1291,1303], 'SkyRight':[1359,1380], 'ObjectAreas':[[1312,1356]]},#lines at 4961 and 5005, no central star visible
##{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0036/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.8e-15], 'OName': 'IPHASXJ214032', 'SkyLeft':[1090,1250], 'SkyRight':[1552,1661], 'ObjectAreas':[[1312,1356]]},#lines at 4961 and 5005, no central star visible
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0037/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-15,0.5e-15], 'OName': 'DSHJ192315', 'SkyLeft':[1185,1214], 'SkyRight':[1687,1748], 'ObjectAreas':[[1277,1299],[1320,1329]]},#lines at 4961 and 5005, pretty compact object
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0038/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0.15e-14], 'OName': 'DSHJ204858', 'SkyLeft':[908,965], 'SkyRight':[1745,1890], 'ObjectAreas':[[1090,1280],[1436,1556]]},#lines at 4860, 4958, 5005, central object gets wider in some wavelength ranges, saturated central? star bleeding into the nebula
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0039/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-15,0.9e-15], 'OName': 'Kn20', 'SkyLeft':[878,964], 'SkyRight':[1570,1606], 'ObjectAreas':[[1296,1332]]},#lines at 6563 and 6582, only on right side of a star (no central star visible?)
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0040/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-15,1.0e-15], 'OName': 'DSHJ195813', 'SkyLeft':[708,799], 'SkyRight':[1603,1655], 'ObjectAreas':[[1214,1257],[1288,1295],[1325,1354],[1383,1394]]},#lines at 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0042/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Kn59', 'SkyLeft':[1090,1185], 'SkyRight':[1412,1511], 'ObjectAreas':[[1285,1303],[1322,1340]]},#compact, lines at 4861, 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0043/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': '', 'SkyLeft':[406,744], 'SkyRight':[1528,1633], 'ObjectAreas':[[1125,1252],[1270,1302],[1337,1503]]},#lines at 4956, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0044/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-2.0e-15,2.0e-15], 'OName': 'Alves1', 'SkyLeft':[286,446], 'SkyRight':[1752,1939], 'ObjectAreas':[[1170,1297],[1324,1471]]},#nothing, without central object
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0044a/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-2.0e-15,2.0e-15], 'OName': 'Alves1', 'SkyLeft':[286,446], 'SkyRight':[1752,1939], 'ObjectAreas':[[1170,1471]]},#nothing, with central object
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0045/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'DSHJ205943', 'SkyLeft':[1090,1143], 'SkyRight':[1348,1403], 'ObjectAreas':[[1291,1304],[1317,1329]]},#compact, lines at 4858, 4956, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0046/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Ou5', 'SkyLeft':[1123,1257], 'SkyRight':[1440,1512], 'ObjectAreas':[[1274,1323],[1342,1385]]},#lines at 4861, 4959, 5005, 5873
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0048/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'DSHJ062355', 'SkyLeft':[811,1028], 'SkyRight':[1677,1790], 'ObjectAreas':[[1113,1135],[1160,1213],[1225,1278],[1301,1332],[1366,1430],[1472,1518],[1524,1561],[1566,1608]]},#lines at 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0049/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ193849', 'SkyLeft':[336,404], 'SkyRight':[1836,1916], 'ObjectAreas':[[944,998],[1012,1038],[1090,1110],[1138,1150],[1168,1220],[1254,1270],[1284,1354],[1382,1412],[1440,1460],[1486,1510],[1542,1552],[1572,1586],[1602,1610],[1626,1652]]},#lines at 6580
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0052/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Ou2', 'SkyLeft':[831,909], 'SkyRight':[1800,1984], 'ObjectAreas':[[1128,1263],[1284,1300],[1322,1327],[1391,1401],[1406,1443],[1448,1483]]},#lines at 4684, 4957, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0055/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.3e-14,0.2e-14], 'OName': 'We2-260', 'SkyLeft':[512,772], 'SkyRight':[1677,1777], 'ObjectAreas':[[871,1038],[1090,1400],[1426,1560],[1567,1657]]},#line at 6582
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0056/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-14,0.6e-14], 'OName': 'Si1-2', 'SkyLeft':[714,954], 'SkyRight':[1524,1578], 'ObjectAreas':[[1167,1292],[1335,1369],[1386,1406],[1419,1436]]},#lines at 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0057/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-15,1.0e-15], 'OName': 'IPHASXJ023538', 'SkyLeft':[719,956], 'SkyRight':[1784,1987], 'ObjectAreas':[[1247,1321],[1370,1449]]},#faint continuum both sides of the central faint object (not included)
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0062/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'DSHJ202907', 'SkyLeft':[926,1034], 'SkyRight':[1600,1847], 'ObjectAreas':[[1227,1274]]},#lines at 4860, 4958, 5005, bright star offset to center (not included)
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0065/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-15,0.1e-15], 'OName': 'Pa30', 'SkyLeft':[839,935], 'SkyRight':[1740,1823], 'ObjectAreas':[[1173,1216],[1323,1344],[1509,1558],[1622,1650],[1697,1737]]},#possible doublet at 6697, 6711 and redshifted at 6737, 6750
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0067/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': '', 'SkyLeft':[394,486], 'SkyRight':[1964,2025], 'ObjectAreas':[[711,766],[796,916],[936,997],[1020,1038],[1131,1208],[1365,1565],[1707,1765],[1799,1822]]},#extended object with lines at 4957, 5005,...
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0068/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-1.0e-15,2.0e-15], 'OName': 'HaWe15', 'SkyLeft':[395,487], 'SkyRight':[1965,2026], 'ObjectAreas':[[712,767],[797,875],[896,917],[938,998],[1021,1038],[1132,1209],[1366,1401],[1422,1470],[1502,1561],[1708,1766],[1800,1823]]},#line at 6562
#####{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0068/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-1.0e-15,2.0e-15], 'OName': 'HaWe15', 'SkyLeft':[1090,1169], 'SkyRight':[1407,1541], 'ObjectAreas':[[1266,1345]]},#line at 6562
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0069/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-2.0e-16,5.0e-15], 'OName': 'We2-5', 'SkyLeft':[710,772], 'SkyRight':[1925,1973], 'ObjectAreas':[[879,1037],[1090,1125],[1143,1233],[1251,1297],[1397,1446],[1514,1711],[1742,1786],[1846,1907]]},#lines at 6546, 6562, 6582
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0070/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-14,0.3e-14], 'OName': 'KLSS2-8', 'SkyLeft':[702,1017], 'SkyRight':[1560,1740], 'ObjectAreas':[[1136,1299],[1325,1409],[1437,1479]]},#line at 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0073/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-14,0.9e-14], 'OName': 'PM1-305', 'SkyLeft':[1149,1211], 'SkyRight':[1625,1663], 'ObjectAreas':[[1283,1303],[1322,1359]]},#lines at 4860, 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0079/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.05e-15,5.0e-15], 'OName': 'Kn33', 'SkyLeft':[778,891], 'SkyRight':[1791,1879], 'ObjectAreas':[[1305,1351],[1376,1418]]},#lines at 4861, 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0080/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.5e-15,1.0e-15], 'OName': 'KKR62', 'SkyLeft':[860,962], 'SkyRight':[1577,1769], 'ObjectAreas':[[1195,1295],[1318,1378]]},#line at 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ191306.1+025248', 'SkyLeft':[429,475], 'SkyRight':[1780,1837], 'ObjectAreas':[[929,972],[992,1038],[1090,1116],[1157,1196],[1234,1288],[1322,1348],[1376,1403],[1471,1526],[1640,1665],[1688,1714]]},#lines at 4861, 5005, 5200
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0003/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.2e-16,1.5e-16], 'OName': 'IPHASXJ190543.8+064413', 'SkyLeft':[717,921], 'SkyRight':[1567,1756], 'ObjectAreas':[[1267,1338]]},#lines at 6550, 6563, 6583, no central object visible
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0004/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.1e-15,0.25e-15], 'OName': 'IPHASXJ191058.9+040350', 'SkyLeft':[840,920], 'SkyRight':[1608,1674], 'ObjectAreas':[[1277,1292],[1309,1317],[1337,1368]]},#lines at 4957, 5005, no central object visible
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0005/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ191104.8+060845', 'SkyLeft':[786,1019], 'SkyRight':[1411,1512], 'ObjectAreas':[[1277,1323],[1345,1362]]},#lines at 4959, 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0007/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-2.5e-16,2.5e-16], 'OName': 'IPHASXJ191707.3+020010', 'SkyLeft':[670,708], 'SkyRight':[1951,2005], 'ObjectAreas':[[1230,1271],[1300,1317],[1351,1385],[1420,1453]]},#lines at 4860 and 5005
{'FName': 'IPHAS_GTC_DATA/DATA4_2018A/GTC11-18AMEX/OB0008/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.02e-14,0.2e-14], 'OName': 'IPHASXJ191918.8+171148', 'SkyLeft':[1090,1186], 'SkyRight':[1451,1557], 'ObjectAreas':[[1251,1284],[1323,1328],[1340,1395]]}#lines at 6548, 6561, 6583
]
i=1
for obs in tab:
    directory = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['FName'][0:obs['FName'].rfind('/')])
    print 'directory = <'+directory+'>'
    print dir(iraf)
#    iraf.user.reduce_gtcmos(indirs=directory)
    image_file = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['FName'])
    try:
        hdulist = pyfits.open(image_file)
    except:
        print 'file <'+image_file+'> not found'
        image_file = image_file[0:image_file.rfind('_x')]+image_file[image_file.rfind('_x')+2:]
        print 'trying <'+image_file+'>'
        try:
            hdulist = pyfits.open(image_file)
        except:
            print 'file <'+image_file+'> not found'

#    print 'type(hdulist) = ',type(hdulist)
    header = hdulist[0].header
    if (obs['OName'] == '') and (header['OBJECT'] != ''):
        obs['OName'] = header['OBJECT']
    try:
        print 'airmass = ',header['AIRMASS']
    except:
        """do nothing"""
    print obs
    wavelength = getWavelength(header)
    image_data = fits.getdata(image_file)
#        print 'type(image_data) = ',type(image_data)
#        print 'type(image_data[1000,1000]) = ',type(image_data[1000,1000])
#        print 'len(image_data) = ',len(image_data)
#        print 'image_data.shape = ',image_data.shape
#        print 'image_data[1551,',obs['ObjectAreas'][0][0]+1,'] = ',image_data[1551,obs['ObjectAreas'][0][0]+1]
    imageMinusSky, imageSky = subtractSky(image_data,obs['SkyLeft'],obs['SkyRight'],3.,2.)
#        print 'imageMinusSky[1551,',obs['ObjectAreas'][0][0]+1,':',obs['ObjectAreas'][0][0]+5,'] = ',imageMinusSky[1551,obs['ObjectAreas'][0][0]+1:obs['ObjectAreas'][0][0]+5]

    obsCols = populateObjectArray(imageMinusSky,obs['ObjectAreas'])
#        print 'obsCols[1551,1:5] = ',obsCols[1551,1:5]
#        print 'obsCols[1000,0] = %.3e' % (obsCols[1000,0])
#        print 'sum(obsCols[1551,:]) = ',np.sum(obsCols[1551,:])

    obsColsSmoothed = boxCarMedianSmooth(obsCols, 0, 5)
#        print 'sum(obsColsSmoothed[1551,:]) = ',np.sum(obsColsSmoothed[1551,:])
#        print 'obsColsSmoothed[1551,1:5] = ',obsColsSmoothed[1551,1:5]
#        print 'imageMinusSky[1551,',obs['ObjectAreas'][0][0]+1,':',obs['ObjectAreas'][0][0]+5,'] = ',imageMinusSky[1551,obs['ObjectAreas'][0][0]+1:obs['ObjectAreas'][0][0]+5]
#        print 'obsCols[1551,1:5] = ',obsCols[1551,1:5]
#        print 'obsColsSmoothed[1551,1:5] = ',obsColsSmoothed[1551,1:5]

    spectrum = np.ndarray(imageSky.shape[0], dtype=np.float32)
    for iRow in range(spectrum.shape[0]):
        spectrum[iRow] = np.sum(obsColsSmoothed[iRow,:])
#            print 'spectrum[',iRow,'] = ',spectrum[iRow]

    plt.plot(wavelength,spectrum)
    plt.xlabel('wavelength [A]')
    yLabel = 'calibrated flux [erg/cm2/s/A]'
    if obs['FName'].find("_cal") < 0:
        print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
        yLabel = 'counts [ADU]'
    plt.ylabel(yLabel)
    if obs['ylim'][1] > 0:
        plt.ylim(obs['ylim'][0],obs['ylim'][1])
    plt.title(obs['OName'])
    plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum.png')
    plt.show()

    # --- mark sky and object areas in original image
    markAreas(image_data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    markAreas(imageMinusSky, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])

    hdulist[0].data = obsCols
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsCols.fits', clobber=True)

    hdulist[0].data = obsColsSmoothed
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsColsSmoothed.fits', clobber=True)

    hdulist[0].data = image_data
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_areasMarked.fits', clobber=True)

    hdulist[0].data = imageSky
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_sky.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    nCols = 0
    for area in obs['ObjectAreas']:
        hdulist[0].data[:,area[0]:area[1]] = obsCols[:,nCols:nCols+area[1]-area[0]]
        nCols += area[1]-area[0]
    markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_not_smoothed.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    nCols = 0
    for area in obs['ObjectAreas']:
        hdulist[0].data[:,area[0]:area[1]] = obsColsSmoothed[:,nCols:nCols+area[1]-area[0]]
        nCols += area[1]-area[0]
    markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_smoothed.fits', clobber=True)

    # --- cut off blue end
    idx = findFirstIdxWithValGT(wavelength, 3700.)
    wavelength = wavelength[idx:]
    spectrum = spectrum[idx:]
#        print 'wavelength.shape = ',wavelength.shape
    print 'spectrum.shape = ',spectrum.shape

    # --- interpolate over 5567-5590 (5577 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 5567.)
    idxH = findFirstIdxWithValGT(wavelength, 5590.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5], 3.0, 3.0)

    # --- interpolate over 5867-5905 (~5885 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 5867.)
    idxH = findFirstIdxWithValGT(wavelength, 5905.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5], 3.0, 3.0)

    # --- interpolate over 6290-6308 (~6299 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 6290.)
    idxH = findFirstIdxWithValGT(wavelength, 6308.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5], 3.0, 3.0)

#    fitsName = image_file[0:image_file.rfind('.')]+'s.fits'
#    write1DFits(spectrum, fitsName)

    # --- continuum normalize spectrum
#    fitsNameOut = fitsName[0:fitsName.rfind('.')]+'c.fits'
#    if os.path.exists(fitsNameOut):
#        os.remove(fitsNameOut)
#    iraf.continuum(fitsName,
#                   fitsNameOut,
#                   ask="yes",
#                   lines="*",
#                   bands="1",
#                   type="ratio",
#                   replace=False,
#                   wavescale=True,
#                   logscale=False,
#                   override=False,
#                   listonly=False,
#                   logfiles="logfile",
#                   interactive=True,
#                   sample="*",
#                   naverage=1,
#                   function="spline3",
#                   order=6,
#                   low_reject=5.,
#                   high_reject=1.,
#                   niterate=10,
#                   grow=1.,
#                   markrej=True)
#    spectrum=read1DFits(fitsNameOut)

#    pymodelfit.fitgui.FitGui(xdata=wavelength, ydata=spectrum, weights=None, model='smoothspline', include_models=None, exclude_models=None, fittype=None, **traits)

    # --- plot cleaned spectrum
    plt.plot(wavelength,spectrum)
    plt.xlabel('wavelength [A]')
    yLabel = 'calibrated flux [erg/cm2/s/A]'
    if obs['FName'].find("_cal") < 0:
        print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
        yLabel = 'counts [ADU]'
    plt.ylabel(yLabel)
    if obs['ylim'][1] > 0:
        plt.ylim(obs['ylim'][0],obs['ylim'][1])
    plt.title(obs['OName'])
    plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum.png')
    plt.show()

    # --- write spectrum
    hdulist[0].data = spectrum
#        hdulist[0].header['NAXIS1'] = header['NAXIS2']
    hdulist[0].header['CRPIX1'] = header['CRPIX2']
    hdulist[0].header['CRVAL1'] = wavelength[0]
    hdulist[0].header['CDELT1'] = header['CDELT2']
    hdulist[0].header['CD1_1'] = header['CD2_2']
    hdulist[0].header['WAT1_001'] = header['WAT2_001']
#        hdulist[0].header[''] = header['']
    obsdate = getDateTime(header['DATE-OBS'])
    datestr = obsdate.strftime('%d%m%y')
    print 'datestr = <'+datestr+'>'
    specOutName = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['OName']+'_GT'+datestr+'.fits')

    hdulist.writeto(specOutName, clobber=True)

    with open(image_file[0:image_file.rfind('.')]+'_spec.dat','w') as f:
        for i in range(len(wavelength)):
            f.write('%.5e %.5e\n' % (wavelength[i], spectrum[i]))
#    STOP
#    plt.imshow(image_data, cmap='gray')
#    plt.colorbar()
#    i+=1
#    if i == 3:
#        STOP
#print tab

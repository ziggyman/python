import paramiko
import os
import math
#from astropy.coordinates import Angle

server = 'tcs74v2.suth.saao.ac.za'
#server = 'shocnhorror.saao'
user='ccd'
passwd='Saaoccd'
remotepath="/home/ccd/live.dat"
#remotepath="/home/ccd/archive/040315.dat"
localpath="live.dat"

def create_ssh_client(server, username, password):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, username=username, password=password)
    return client

def read_remote_file(server, user, password, remotepath):
    print("Creating ssh client to {}@{}".format(user, server))
    try:
        ssh = create_ssh_client(server, user, password)
        sftp = ssh.open_sftp()
        print("Reading file {}:{}".format(server, remotepath))
        with sftp.open(remotepath) as remote_file:
            lines = remote_file.read()
        sftp.close()
        ssh.close()
    except (paramiko.ssh_exception.BadHostKeyException,
            paramiko.ssh_exception.AuthenticationException,
            paramiko.ssh_exception.BadAuthenticationType,
            paramiko.ssh_exception.ChannelException, 
            paramiko.ssh_exception.SSHException) as e:
        raise Exception("Paramiko exception {}".format(e))
    except IOError as e:
        raise Exception("Error reading file {}:{}: {}".format(server, remotepath, str(e)))
    return lines

def write_remote_file(server, user, password, remotepath, lines):
    print("Creating ssh client to {}@{}".format(user, server))
    try:
        ssh = create_ssh_client(server, user, password)
        sftp = ssh.open_sftp()
        print("Writing file {}:{}".format(server, remotepath))
        remote_file=sftp.file(remotepath, "w", -1)        
#        with sftp.open(remotepath) as remote_file:
        remote_file.write(lines)
        remote_file.flush()
        sftp.close()
        ssh.close()
    except (paramiko.ssh_exception.BadHostKeyException,
            paramiko.ssh_exception.AuthenticationException,
            paramiko.ssh_exception.BadAuthenticationType,
            paramiko.ssh_exception.ChannelException, 
            paramiko.ssh_exception.SSHException) as e:
        raise Exception("Paramiko exception {}".format(e))
    except IOError as e:
        raise Exception("Error reading file {}:{}: {}".format(server, remotepath, str(e)))
#    return lines



"""Get the last line."""
def get_last_line(lines):
    for line in lines.split("\n"):
        if len(line.split(" ")) > 1:
            lastline = line
    return lastline.strip()

def split_line(line):
    tokens = line.split()
    tcs = {}
    tcs["JULIAN_DATE"] = tokens[0]
    tcs["TELRA"] = ":".join(tokens[1:4])
    tcs["TELDEC"] = ":".join(tokens[4:7])
    tcs["SECZ"] = tokens[7]
    tcs["ST"] = ":".join(tokens[8:11])
    tcs["TELFOCUS"] = tokens[11]
    tcs["INSTANGL"] = tokens[12]
    tcs["DOMEPOS"] = tokens[13]
    tcs["DOMEREQD"] = tokens[14]
    tcs["TARGRA"] = ":".join(tokens[15:18])
    tcs["TARGDEC"] = ":".join(tokens[18:21])
    tcs["RAZERO"] = tokens[21]
    tcs["DECZERO"] = tokens[22]
    tcs["RAWRA"] = ":".join(tokens[24:27])
    tcs["RAWDEC"] = ":".join(tokens[27:30])
    return tcs

def get_tcs_info(line):
    tcs = split_line(line)
    tcs["AIRMASS"] = tcs["SECZ"]
    tcs["ZD"] = str(math.acos(1/float(tcs["SECZ"])))
    st_hours = "{} {} {} hours".format(*tcs["ST"].split(":"))
    ra_hours = "{} {} {} hours".format(*tcs["TELRA"].split(":"))
    st = Angle(st_hours)
    ra = Angle(ra_hours)
    ha = (st - ra).to_string(sep=":")
    tcs["HA"] = ha
    return tcs

def show_tcs_info(tcs):
    fmtstr = """
Telescope Control System information:

Julian date                = {JULIAN_DATE}
Right ascension            = {TELRA}
Declination                = {TELDEC}
secZ                       = {SECZ}
Sidereal time              = {ST}
Telescope focus            = {TELFOCUS}
Instrument angle           = {INSTANGL}
Dome position              = {DOMEPOS}
Required dome position     = {DOMEREQD}
Target right ascension     = {TARGRA}
Target declination         = {TARGDEC}
Dome position              = {DOMEPOS}
Right ascension zero       = {RAZERO}
Declination zero           = {DECZERO}
Right ascension (raw)      = {RAWRA}
Declination (raw)          = {RAWDEC}
Airmass                    = {AIRMASS}
Zenith distance            = {ZD}
Hour angle                 = {HA}
"""
    print(fmtstr.format(**tcs))

if __name__ == '__main__':
    # fetch_file(server, user, passwd, remotepath, localpath)
    
    # with open(localpath, "r") as fin:
    #     lastline = parse_file(fin)
    try:
        lines = read_remote_file(server, user, passwd, remotepath)
    except Exception as e:
        print("Caught exception: {}".format(e))
    lastline = get_last_line(lines)
    tcs = get_tcs_info(lastline)
    show_tcs_info(tcs)

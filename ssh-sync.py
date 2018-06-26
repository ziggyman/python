#! /usr/bin/env python
import errno
import numpy as np
import os
from subprocess import Popen    # now we can reference Popen
import sys

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

class listContent:
    prefixMac = ''
    prefixDesktop = ''
    ipDesktop = ''

    def isLineEmpty(self, line):
        return len(line.strip()) == 0

    def readContentFile(self, fname):
        text_file = open(fname, "r")
        lines = text_file.readlines()

#        """remove empty lines"""
#        lines = list(filter(len, lines))

        linesOut = []
        dir = None
        lineOut = None
        for i in range(len(lines)):
            if self.isLineEmpty(lines[i]):
                dir = None
            else:
                lines[i] = lines[i].split(' ')

                """remove empty strings from line"""
                lines[i] = list(filter(len, lines[i]))
                print 'lines[',i,'] = ',lines[i]

                lastCol = len(lines[i])-1
                print 'lastCol = ',lastCol

                """remove \n from end of last element in lines[i]"""
                lines[i][lastCol] = lines[i][lastCol][0:len(lines[i][lastCol])-1]
    #            print "lines[',i,']='",lines[i],"' != ''"
    #            print "len(lines[',i,'])='",len(lines[i])

                lastCharPos = len(lines[i][lastCol])-1
                print 'lastCharPos = ',lastCharPos
                if lastCharPos < 0:
                    lastCol -= 1
                    lastCharPos = len(lines[i][lastCol])-1
    #            print 'lines[i][lastCol] = ',lines[i][lastCol]
    #            print 'new lines[',i,'] = ',lines[i]

                if lines[i][lastCol][lastCharPos] == ":":
                    dir = lines[i][0]
                    for dirCol in np.arange(1,len(lines[i])):
                        dir += ' '+lines[i][dirCol]
                    dir = dir[0:len(dir)-1]
                    dir += '/'
                    print 'dir = ',type(dir),': <',dir,'>'
                    dir = dir.replace(self.prefixDesktop+'/','')
                    dir = dir.replace(self.prefixMac+'/','')
                    dir = dir.replace('(','\(')
                    dir = dir.replace(')','\)')
                    dir = dir.replace(' ','\ ')
#                    print 'dir = <',dir,'>'
    #                print 'new dir = ',dir
                elif (dir is not None) and (len(lines[i]) >= 9) and (lines[i][0][0] != 'd'):
    #                print 'dir = <',dir,'>, fileName = <',lines[i][lastCol],'>, size = int(',lines[i][4],')'
                    fileName = lines[i][8]
                    print 'fileName = <'+fileName+'>'
                    for iCol in np.arange(9, lastCol+1, 1):
                        fileName += ' '+lines[i][iCol]
                    fileName = fileName.replace('(','\(')
                    fileName = fileName.replace(')','\)')
                    fileName = fileName.replace(' ','\ ')
                    print 'fileName = <'+fileName+'>'
#                    if ('xy' not in fileName):
                    lineOut = {'dir': dir, 'fileName':fileName, 'size':int(lines[i][4])}
                    print 'lineOut = ',lineOut
                    linesOut.append(lineOut)
        return linesOut

    def findFile(self, dir, fileName, fileDicts):
        dirListMac = list(lineMac['dir'] for lineMac in fileDicts)
#        print 'dirListMac = ',dirListMac

        for fileDict in fileDicts:
            if dir == fileDict['dir']:
#                print 'dir <',dir,'> found in linesMac'
                if fileName == fileDict['fileName']:
                    print 'file <'+dir+'/'+fileName+'> found'
                    return fileDict
        return None

def main(argv):
    ipDesktop = '192.168.43.31'

#    dirMac = "/Volumes/yoda/azuri/data/galaxia/ubv_Vlt21.5_1.0.bak"
    dirMac = "/Volumes/external/azuri/Photos"
    contentFileMac = dirMac+".content"

#    dirDesktop = "/Volumes/yoda/azuri/data/galaxia/ubv_Vlt21.5_1.0.bak"
    dirDesktop = "/home/azuri/daten/bilder"
#    contentFileDesktop = dirDesktop+".content.desktop"
    contentFileDesktop = dirMac+".content.desktop"

    command = 'ls -lR '+dirMac+' > '+contentFileMac
    if os.system(command):
        print 'command <'+command+'> failed'
        STOP

    command = 'ssh azuri@'+ipDesktop+' ls -lR '+dirDesktop+' > '+contentFileDesktop
    if os.system(command):
        print 'command <'+command+'> failed'
        STOP

    a = listContent()
    a.prefixMac = dirMac
    a.prefixDesktop = dirDesktop
    a.ipDesktop = ipDesktop

    linesMac = a.readContentFile(contentFileMac)
#    print 'linesMac = ',linesMac
    print 'linesMac[1] = ',linesMac[1]
#    STOP

    linesDesktop = a.readContentFile(contentFileDesktop)
#    STOP

    with open('copyErrors', 'w+') as f:

        for line in linesDesktop:
            fileIsGood = False
            macLine = a.findFile(line['dir'], line['fileName'], linesMac)
            if macLine != None:
                if macLine['size'] == line['size']:
                    fileIsGood = True
                    print 'file <'+line['dir']+line['fileName']+'> is good'
                else:
                    print 'file <'+line['dir']+line['fileName']+'> is NOT good'
    #                STOP
            if not fileIsGood:
                print 'copying <',line['dir']+'/'+line['fileName'],'>'
                mkdir_p(a.prefixMac+'/'+line['dir'])
                command = 'scp azuri@'+a.ipDesktop+':'+a.prefixDesktop+'/'+line['dir'].replace('\\', '\\\\\\')+line['fileName'].replace('\\', '\\\\\\')+' '+a.prefixMac+'/'+line['dir']#+line['fileName']
    #            print 'command = <',command,'>'
                #process = Popen([command,'/c '+commandArguments])
                #process.wait()
                if os.system(command):
                    print 'command <'+command+'> failed'
    #                STOP
                    f.write(line['dir']+'/'+line['fileName']+'\n')

if __name__ == '__main__':
    main(sys.argv)

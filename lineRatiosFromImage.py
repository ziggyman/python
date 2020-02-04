import cv2
import fitz
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os
import PyPDF2
from PIL import Image

from drUtils import subtractBackground

def doQuentinsImages():
    imA = '/Users/azuri/daten/uni/HKU/line_ratios/Wyse-1942-plateXII.png'
    imB = '/Users/azuri/daten/uni/HKU/line_ratios/Wyse-1942-plateXIII.png'

    im = Image.open(imA) # Can be many different formats.
    pix = im.load()
    print(im.size)  # Get the width and hight of the image for iterating over
    print(pix[0,0])  # Get the RGBA Value of the a pixel of an image
    plt.imshow(im)
    plt.title('Plate XII')
    plt.xlabel('column')
    plt.ylabel('row')
    plt.show()

    print('pix[',im.size[0]-1,',',im.size[1]-1,'] = ',pix[im.size[0]-1,im.size[1]-1])

    ratio5007Hbeta = []
    ratio50074959 = []

    indicesToIgnore = []
    for i in np.arange(50,66,1):
        indicesToIgnore.append(i)
    for i in np.arange(88,102,1):
        indicesToIgnore.append(i)
    for i in np.arange(106,122,1):
        indicesToIgnore.append(i)

    for row in np.arange(20,75,1):
        r = []
        for col in np.arange(0,im.size[0],1):
    #        print('col = ',col)
            r.append(pix[int(col),int(row)][0])
        r = 255-np.array(r)
        rmb = subtractBackground(np.arange(r.shape[0]),r,2,indicesToIgnore)
        if row < 30:#> 64:
            plt.plot(rmb, label='row '+str(row))
        if (np.amax(r[106:122]) < 254) and (np.amax(r[50:66]) < 254):
            ratio5007Hbeta.append(np.amax(rmb[106:122]) / np.amax(rmb[50:66]))
        else:
            ratio5007Hbeta.append(0)
        if (np.amax(r[106:122]) < 254) and (np.amax(r[88:102]) < 254):
            ratio50074959.append(np.amax(rmb[106:122]) / np.amax(rmb[88:102]))
        else:
            ratio50074959.append(0)
        print('row = ',row,': 5007 = ',rmb[106:122])
        print('row = ',row,': Hbeta = ',rmb[50:66])
        print('row = ',row,': 4959 = ',rmb[88:102])
    plt.xlabel('column')
    plt.ylabel('counts')
    plt.title('Plate XII')
    plt.legend()
    plt.show()
    plt.plot(ratio5007Hbeta,label='plateXII 5007/Hbeta')
    #plt.legend()
    plt.xlabel('row - 20')
    plt.ylabel('5007/Hbeta')
    plt.title('Plate XII')
    plt.show()
    plt.plot(ratio50074959,label='plateXII 5007/4959')
    #plt.legend()
    plt.xlabel('row - 20')
    plt.ylabel('5007/4959')
    plt.title('Plate XII')
    plt.show()

    #img=mpimg.imread(imB)
    #print('img = ',img)
    #imgplot = plt.imshow(img)
    #plt.show()
    im = Image.open(imB) # Can be many different formats.
    img = im.load()
    imgplot = plt.imshow(im)
    plt.title('Plate XIII')
    plt.xlabel('column')
    plt.ylabel('row')
    plt.show()

    ratio5007Hbeta = []
    ratio50074959 = []

    indicesToIgnore = []
    for i in np.arange(10,33,1):
        indicesToIgnore.append(i)
    for i in np.arange(128,145,1):
        indicesToIgnore.append(i)
    for i in np.arange(177,207,1):
        indicesToIgnore.append(i)

    for row in np.arange(15,79,1):
        r = []
        for col in np.arange(223,453,1):
            r.append(img[int(col),int(row)][0])
        r = 255 - np.array(r)
        rmb = subtractBackground(np.arange(r.shape[0]),r,2,indicesToIgnore)

        #if row < 25:#> 68:
        #    plt.plot(rmb, label='row '+str(row))
        plt.plot(rmb)
        if (np.amax(r[177:207]) < 254) and (np.amax(r[10:33]) < 254):
            ratio5007Hbeta.append(np.amax(rmb[177:207]) / np.amax(rmb[10:33]))
        else:
            ratio5007Hbeta.append(0)
        if (np.amax(r[177:207]) < 254) and (np.amax(r[128:145]) < 254):
            ratio50074959.append(np.amax(rmb[177:207]) / np.amax(rmb[128:145]))
        else:
            ratio50074959.append(0)
        print('row = ',row,': 5007 = ',rmb[177:207])
        print('row = ',row,': Hbeta = ',rmb[10:33])
        print('row = ',row,': 4959 = ',rmb[128:145])
    plt.xlabel('column - 223')
    plt.ylabel('counts')
    plt.title('Plate XIII')
    plt.legend()
    plt.show()

    plt.plot(ratio5007Hbeta,label='plateXIII 5007/Hbeta')
    plt.xlabel('row - 15')
    plt.ylabel('5007/Hbeta')
    plt.title('Plate XIII')
    #plt.legend()
    plt.show()
    plt.plot(ratio50074959,label='plateXIII 5007/4959')
    #plt.legend()
    plt.xlabel('row - 15')
    plt.ylabel('5007/4959')
    plt.title('Plate XIII')
    plt.show()

def extractImagesFromPdf(pdfFileName):
    path = pdfFileName[0:pdfFileName.rfind('/')]
    print('path = <'+path+'>')
    if False:
        doc = fitz.open(pdfFileName)
        for i in range(len(doc)):
            for img in doc.getPageImageList(i):
                xref = img[0]
                pix = fitz.Pixmap(doc, xref)
                if pix.n < 5:       # this is GRAY or RGB
                    pix.writePNG("p%s-%s.png" % (i, xref))
                else:               # CMYK: convert to RGB first
                    pix1 = fitz.Pixmap(fitz.csRGB, pix)
                    pix1.writePNG("p%s-%s.png" % (i, xref))
                    pix1 = None
                pix = None

    input1 = PyPDF2.PdfFileReader(open(pdfFileName, "rb"))
    outFiles = []
    for iPage in [0,1]:
        page0 = input1.getPage(iPage)
        xObject = page0['/Resources']['/XObject'].getObject()

        for obj in xObject:
            if xObject[obj]['/Subtype'] == '/Image':
                size = (xObject[obj]['/Width'], xObject[obj]['/Height'])
                data = xObject[obj].getData()
                if xObject[obj]['/ColorSpace'] == '/DeviceRGB':
                    mode = "RGB"
                else:
                    mode = "P"

                outFile = ''
                if xObject[obj]['/Filter'] == '/FlateDecode':
                    img = Image.frombytes(mode, size, data)
                    outFile = os.path.join(path, obj[1:] + ".png")
                    img.save(outFile)
                elif xObject[obj]['/Filter'] == '/DCTDecode':
                    outFile = os.path.join(path, obj[1:] + ".jpg")
                    img = open(outFile, "wb")
                    img.write(data)
                    img.close()
                elif xObject[obj]['/Filter'] == '/JPXDecode':
                    outFile = os.path.join(path, obj[1:] + ".jp2")
                    img = open(outFile, "wb")
                    img.write(data)
                    img.close()
                outFiles.append(outFile)
    return outFiles

def findPixGT(pixels, threshold):
    out = np.where(np.array(pixels) > threshold)
#    print(out)
    return out

def doFromPdf():
    pdfFileName = "/Users/azuri/daten/uni/HKU/line_ratios/ApJ95PLATExii001.pdf"
    imFiles = extractImagesFromPdf(pdfFileName)
    print('imFiles = ',imFiles)
    imXRange = None
    imYRange = None
    ignoreColsA = None
    ignoreColsB = None
    ignoreColsC = None
    for iIm in range(len(imFiles)):
        if iIm == 0:
            imXRange = [4325,4650]
            imYRange = [2809,2990]
            ignoreColsA = np.arange(4386,4439)
            ignoreColsB = np.arange(4512,4556)
            ignoreColsC = np.arange(4567,4615)
        else:
            imXRange = [1315,2032]
            imYRange = [3207,3442]
            ignoreColsA = np.arange(1356,1406,1)
            ignoreColsB = np.arange(1723,1770,1)
            ignoreColsC = np.arange(1886,1944,1)
        imName = imFiles[iIm]
        im = Image.open(imName) # Can be many different formats.
        pix = im.load()
        print(im.size)  # Get the width and hight of the image for iterating over
        print('pix[4429,2783] = ',pix[4429,2783])  # Get the RGBA Value of the a pixel of an image
        print('pix[4412,2866] = ',pix[4412,2866])  # Get the RGBA Value of the a pixel of an image

        plt.imshow(im)
        plt.title('Plate XII')
        plt.xlabel('column')
        plt.ylabel('row')
        plt.show()

        ratio5007Hbeta = []
        ratio50074959 = []

        indicesToIgnore = []
        for i in (ignoreColsA-imXRange[0]):
            indicesToIgnore.append(i)
        for i in (ignoreColsB-imXRange[0]):
            indicesToIgnore.append(i)
        for i in (ignoreColsC-imXRange[0]):
            indicesToIgnore.append(i)
        print('indicesToIgnore = ',indicesToIgnore)

        greyPix = np.ndarray(shape=(imXRange[1]-imXRange[0], imYRange[1]-imYRange[0]))
        for row in np.arange(imYRange[0],imYRange[1],1):
            for col in np.arange(imXRange[0],imXRange[1],1):
                greyPix[col-imXRange[0], row-imYRange[0]] = np.sum(pix[int(col),int(row)])
        plt.imshow(greyPix)
        plt.show()

        greyPixFiltered = cv2.blur(greyPix,(5,5))
        plt.imshow(greyPixFiltered)
        plt.show()
        greyPixFiltered = cv2.blur(greyPixFiltered,(5,5))
        plt.imshow(greyPixFiltered)
        plt.show()

        pixelsAboveThreshold5007 = []
        pixelsAboveThreshold4959 = []
        pixelsAboveThresholdHbeta = []
        for row in np.arange(0,imYRange[1]-imYRange[0],1):
            r = []
            for col in np.arange(0,imXRange[1]-imXRange[0],1):
        #        print('col = ',col)
                print('pix[',col,',',row,'] = ',pix[int(col),int(row)],': sum = ',np.sum(pix[int(col),int(row)]))
                r.append(greyPixFiltered[int(col),int(row)])
            print('max(r) = ',np.amax(r))
#            plt.plot(r)
#            plt.show()
            r = (3*255)-np.array(r)
            #r = np.array(r)
            rmb = subtractBackground(np.arange(r.shape[0]),r,3,indicesToIgnore)
            pixelsAboveThreshold5007.append(findPixGT(rmb[ignoreColsC-imXRange[0]], 50.))
            pixelsAboveThresholdHbeta.append(findPixGT(rmb[ignoreColsA-imXRange[0]], 50.))
            pixelsAboveThreshold4959.append(findPixGT(rmb[ignoreColsB-imXRange[0]], 50.))
            print('rmb.size = ',rmb.size)
            plt.plot(np.arange(0,rmb.size,1)+imXRange[0],rmb, label='row '+str(row+imYRange[0]))
            if (row % 10) == 9:
                plt.xlabel('column')
                plt.ylabel('counts')
                plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
                plt.legend()
                plt.show()
            if (np.amax(r[ignoreColsC-imXRange[0]]) < 3*254) and (np.amax(r[ignoreColsA-imXRange[0]]) < 3*254):
                ratio5007Hbeta.append(np.amax(rmb[ignoreColsC-imXRange[0]]) / np.amax(rmb[ignoreColsA-imXRange[0]]))
            else:
                ratio5007Hbeta.append(0)
            if (np.amax(r[ignoreColsC-imXRange[0]]) < 3*254) and (np.amax(r[ignoreColsB-imXRange[0]]) < 3*254):
                ratio50074959.append(np.amax(rmb[ignoreColsC-imXRange[0]]) / np.amax(rmb[ignoreColsB-imXRange[0]]))
            else:
                ratio50074959.append(0)
            print('row = ',row,': 5007 = ',rmb[ignoreColsC-imXRange[0]])
            print('row = ',row,': Hbeta = ',rmb[ignoreColsA-imXRange[0]])
            print('row = ',row,': 4959 = ',rmb[ignoreColsB-imXRange[0]])
#            STOP

        plt.xlabel('column')
        plt.ylabel('counts')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plt.legend()
        plt.show()

        plt.plot(np.arange(0,len(ratio5007Hbeta),1)+imYRange[0],ratio5007Hbeta,label='5007/Hbeta')
        #plt.legend()
        plt.xlabel('row')
        plt.ylabel('5007/Hbeta')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plt.show()

        plt.plot(np.arange(0,len(ratio50074959),1)+imYRange[0],ratio50074959,label='5007/4959')
        #plt.legend()
        plt.xlabel('row')
        plt.ylabel('5007/4959')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plt.show()

        ratio5007HbetaNPix = []
        ratio50074959NPix = []
        for row in np.arange(0,imYRange[1]-imYRange[0],1):
            print('pixelsAboveThreshold5007[',row,'] = ',len(pixelsAboveThreshold5007[row]),': ',pixelsAboveThreshold5007[row])
            print('pixelsAboveThreshold5007[',row,'][0] = ',len(pixelsAboveThreshold5007[row][0]),': ',pixelsAboveThreshold5007[row][0])
            print('pixelsAboveThreshold4959[',row,'] = ',len(pixelsAboveThreshold4959[row]),': ',pixelsAboveThreshold4959[row])
            print('pixelsAboveThresholdHbeta[',row,'] = ',len(pixelsAboveThresholdHbeta[row]),': ',pixelsAboveThresholdHbeta[row])
            ratio5007HbetaNPix.append(len(pixelsAboveThreshold5007[row][0]) / len(pixelsAboveThresholdHbeta[row][0]) if len(pixelsAboveThresholdHbeta[row][0])>0 else 0)
            ratio50074959NPix.append(len(pixelsAboveThreshold5007[row][0]) / len(pixelsAboveThreshold4959[row][0]) if len(pixelsAboveThreshold4959[row][0])>0 else 0)

        plt.plot(np.arange(imYRange[0],imYRange[1],1),ratio5007HbetaNPix)
        plt.xlabel('row')
        plt.ylabel('nPix(5007)/nPix(Hbeta)')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plt.show()

        plt.plot(np.arange(imYRange[0],imYRange[1],1),ratio50074959NPix)
        plt.xlabel('row')
        plt.ylabel('nPix(5007)/nPix(4959)')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plt.show()

def main():
    doFromPdf()

if __name__ == '__main__':
    main()

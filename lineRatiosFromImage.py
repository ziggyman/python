import cv2
#import fitz
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os
import PyPDF2
from PIL import Image
from scipy.optimize import curve_fit

from drUtils import subtractBackground,gauss
import mympfit# import MPFitTwoGaussLim

nPages = 1#2
negative = False

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

#    print('pix[',im.size[0]-1,',',im.size[1]-1,'] = ',pix[im.size[0]-1,im.size[1]-1])

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
#        print('row = ',row,': 5007 = ',rmb[106:122])
#        print('row = ',row,': Hbeta = ',rmb[50:66])
#        print('row = ',row,': 4959 = ',rmb[88:102])
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
#        print('row = ',row,': 5007 = ',rmb[177:207])
#        print('row = ',row,': Hbeta = ',rmb[10:33])
#        print('row = ',row,': 4959 = ',rmb[128:145])
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
#    if False:
#        doc = fitz.open(pdfFileName)
#        for i in range(len(doc)):
#            for img in doc.getPageImageList(i):
#                xref = img[0]
#                pix = fitz.Pixmap(doc, xref)
#                if pix.n < 5:       # this is GRAY or RGB
#                    pix.writePNG("p%s-%s.png" % (i, xref))
#                else:               # CMYK: convert to RGB first
#                    pix1 = fitz.Pixmap(fitz.csRGB, pix)
#                    pix1.writePNG("p%s-%s.png" % (i, xref))
#                    pix1 = None
#                pix = None

    input1 = PyPDF2.PdfFileReader(open(pdfFileName, "rb"))
    outFiles = []
    for iPage in range(nPages):
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

def interpolate(image, rowStart, rowEndOffset):
#    print('interpolate: image.shape = ',image.shape)
    colStart = 0# if rowEndOffset >= 0 else 0-int(rowEndOffset)+1
    colEnd = image.shape[0]#-(int(rowEndOffset) if rowEndOffset >= 0 else 0) - 1
#    print('interpolate: rowStart = ',rowStart,', rowEndOffset = ',rowEndOffset)
#    print('interpolate: colStart = ',colStart,', colEnd = ',colEnd)
    rowOut = np.zeros(colEnd - colStart)
    for col in np.arange(colStart,colEnd,1):
        offset = rowEndOffset * col / image.shape[0]
#        print('col = ',col,': offset = ',offset)
#        print('col = ',col,': image[col, rowStart + int(offset)] = ',image[col, rowStart + int(offset)])
#        print('col = ',col,': image[col, rowStart + int(offset) + ',(1 if rowEndOffset >= 0 else -1),'] = ',image[col, rowStart + int(offset) + (1 if rowEndOffset >= 0 else -1)])
        rowOut[col] = image[col, rowStart + int(offset)] + (
                     (image[col, rowStart + int(offset) + (1 if rowEndOffset >= 0 else -1)] - image[col, rowStart + int(offset)])
                     * (offset - int(offset)))
#        print('col = ',col,': rowOut[col] = ',rowOut[col])
    return rowOut

def doFromPdf():
    pdfFileName = "/Users/azuri/daten/uni/HKU/line_ratios/ApJ140Aler001.pdf"
#    pdfFileName = "/Users/azuri/daten/uni/HKU/line_ratios/ApJ95PLATExii001.pdf"
    path = pdfFileName[0:pdfFileName.rfind('/')]
    imFiles = extractImagesFromPdf(pdfFileName)
    print('imFiles = ',len(imFiles),': ',imFiles)
    imXRange = None
    imYRange = None
    ignoreColsA = None
    ignoreColsB = None
    ignoreColsC = None
    rowOffsetHbetaO5007 = None
    for iIm in range(len(imFiles)):
        if iIm == 0:
            imXRange = [3045,4244]
            imYRange = [1100,1192]
            ignoreColsA = np.arange(3162,3201)
            ignoreColsB = np.arange(3842,3877)
            ignoreColsC = np.arange(4172,4212)
            rowOffsetHbetaO5007 = -17
#            imXRange = [4325,4650]
#            imYRange = [2809,2990]
#            ignoreColsA = np.arange(4386,4439)
#            ignoreColsB = np.arange(4512,4556)
#            ignoreColsC = np.arange(4567,4615)
#            rowOffsetHbetaO5007 = 0
        else:
            imXRange = [1315,2032]
            imYRange = [3207,3442]
            ignoreColsA = np.arange(1356,1406,1)
            ignoreColsB = np.arange(1723,1770,1)
            ignoreColsC = np.arange(1886,1944,1)
            rowOffsetHbetaO5007 = 3
        imName = imFiles[iIm]
        im = Image.open(imName) # Can be many different formats.
        pix = im.load()
        print(im.size)  # Get the width and hight of the image for iterating over
        print('pix[4192,1139] = ',pix[4192,1139])  # Get the RGBA Value of the a pixel of an image
        print('pix[4092,1139] = ',pix[4092,1139])  # Get the RGBA Value of the a pixel of an image

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
        plt.imshow(np.transpose(greyPix), cmap='gray', vmin=0, vmax=255*3)
        plt.colorbar()
        plt.xlabel('column - %d' % (imXRange[0]))
        plt.ylabel('row - %d' % (imYRange[0]))
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_grey.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()

        greyPixFiltered = cv2.blur(greyPix,(5,5))
        plt.imshow(np.transpose(greyPixFiltered), cmap='gray', vmin=0, vmax=255*3)
        plt.colorbar()
        plt.xlabel('column - %d' % (imXRange[0]))
        plt.ylabel('row - %d' % (imYRange[0]))
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_smoothed.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()

        greyPixFiltered = cv2.blur(greyPixFiltered,(5,5))
        plt.imshow(np.transpose(greyPixFiltered), cmap='gray', vmin=0, vmax=255*3)
        plt.colorbar()
        plt.xlabel('column - %d' % (imXRange[0]))
        plt.ylabel('row - %d' % (imYRange[0]))
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_smoothed_twice.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()

        pixelsAboveThreshold5007 = []
        pixelsAboveThreshold4959 = []
        pixelsAboveThresholdHbeta = []
        rowStart = 0 if rowOffsetHbetaO5007 >= 0 else 0-rowOffsetHbetaO5007
        rowEnd = imYRange[1]-imYRange[0]-rowOffsetHbetaO5007-1 if rowOffsetHbetaO5007 >= 0 else imYRange[1]-imYRange[0]-1
        for row in np.arange(rowStart,rowEnd,1):
            r = interpolate(greyPixFiltered[0:imXRange[1]-imXRange[0],0:imYRange[1]-imYRange[0]], row, rowOffsetHbetaO5007)
#            for col in np.arange(0,imXRange[1]-imXRange[0],1):
        #        print('col = ',col)
#                print('pix[',col,',',row,'] = ',pix[int(col),int(row)],': sum = ',np.sum(pix[int(col),int(row)]))
#                r.append(greyPixFiltered[int(col),int(row)])
#                rplus1.append(greyPixFiltered[int(col),int(row+1)])
#                rplus2.append(greyPixFiltered[int(col),int(row+2)])
#                rplus3.append(greyPixFiltered[int(col),int(row+3)])
            print('max(r) = ',np.amax(r))
#            plt.plot(r)
#            plt.show()
            if negative:
                r = (3*255)-np.array(r)
            indices = np.arange(r.shape[0])
            print('indices = ',indices)
            rmb = subtractBackground(indices,r,4,indicesToIgnore)

            xiA = None
            xiB = None
            xiC = None
            sigma = 5.25
            yFitA = None
            yFitB = None
            yFitC = None
            try:
                xiA = indices[ignoreColsA-imXRange[0]]
                print('xiA = ',xiA.shape,': ',xiA)
                yi = rmb[xiA]
                print('yi = ',yi.shape,': ',yi)
                x0=xiA[int(xiA.shape[0]/2.)]
                print('x0 = ',x0)
                xiAToFit = []
                yiToFit = []
                for i in np.arange(0,len(xiA),1):
                    if xiA[i] < x0-(sigma):
                        xiAToFit.append(xiA[i])
                        yiToFit.append(yi[i])
                    if xiA[i] > x0+(sigma):
                        xiAToFit.append(xiA[i])
                        yiToFit.append(yi[i])
                popt,pcov = curve_fit(gauss,
                                      xiAToFit,
                                      yiToFit,
                                      p0=[np.amax(yi),
                                          x0,
                                          sigma,
                                          0.
                                         ]
                                     )
                print('popt = ',popt)
                yFitA = gauss(xiA,*popt)
            except Exception as e:
                print(e)
                continue

            try:
                xiB = indices[ignoreColsB-imXRange[0]]
                print('xiB = ',xiB.shape,': ',xiB)
                yi = rmb[xiB]
                print('yi = ',yi.shape,': ',yi)
                x0=xiB[int(xiB.shape[0]/2.)]
                print('x0 = ',x0)
                xiBToFit = []
                yiToFit = []
                for i in np.arange(0,len(xiB),1):
                    if xiB[i] < x0-(sigma):
                        xiBToFit.append(xiB[i])
                        yiToFit.append(yi[i])
                    if xiB[i] > x0+(sigma):
                        xiBToFit.append(xiB[i])
                        yiToFit.append(yi[i])
                popt,pcov = curve_fit(gauss,
                                      xiBToFit,
                                      yiToFit,
                                      p0=[np.amax(yi),
                                          x0,
                                          sigma,
                                          0.
                                         ]
                                     )
                print('popt = ',popt)
                yFitB = gauss(xiB,*popt)
            except Exception as e:
                print(e)
                continue

            try:
                xiC = indices[ignoreColsC-imXRange[0]]
                print('xiC = ',xiC.shape,': ',xiC)
                yi = rmb[xiC]
                print('yi = ',yi.shape,': ',yi)
                x0=xiC[int(xiC.shape[0]/2.)]
                print('x0 = ',x0)
                xiCToFit = []
                yiToFit = []
                for i in np.arange(0,len(xiC),1):
                    if xiC[i] < x0-(sigma):
                        xiCToFit.append(xiC[i])
                        yiToFit.append(yi[i])
                    if xiC[i] > x0+(sigma):
                        xiCToFit.append(xiC[i])
                        yiToFit.append(yi[i])
                popt,pcov = curve_fit(gauss,
                                      xiCToFit,
                                      yiToFit,
                                      p0=[np.amax(yi),
                                          x0,
                                          sigma,
                                          0.
                                         ]
                                     )
                print('popt = ',popt)
                yFitC = gauss(xiC,*popt)
            except Exception as e:
                print(e)
                continue


            pixelsAboveThreshold5007.append(findPixGT(rmb[ignoreColsC-imXRange[0]], 50.))
            pixelsAboveThresholdHbeta.append(findPixGT(rmb[ignoreColsA-imXRange[0]], 50.))
            pixelsAboveThreshold4959.append(findPixGT(rmb[ignoreColsB-imXRange[0]], 50.))
            print('rmb.size = ',rmb.size)
            plt.plot(np.arange(0,rmb.size,1)+imXRange[0],rmb, label='row '+str(row+imYRange[0]))
            if yFitA is not None:
                plt.plot(xiA+imXRange[0],yFitA)
            if yFitB is not None:
                plt.plot(xiB+imXRange[0],yFitB)
            if yFitC is not None:
                plt.plot(xiC+imXRange[0],yFitC)
            if ((row % 10) == 9) or (row == rowEnd-1):
                plt.xlabel('column')
                plt.ylabel('counts')
                plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
                if iIm == 0:
                    plt.legend(loc='lower left', bbox_to_anchor=(0.35, 0.25))
                else:
                    plt.legend(loc='lower left', bbox_to_anchor=(0.2, 0.33))
                plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
                plotName += '_row%d-%d.png' % (imYRange[0]+row-9,imYRange[0]+row)
                plt.savefig(os.path.join(path,plotName))
                plt.show()
#                STOP
            if (np.amax(r[ignoreColsC-imXRange[0]]) < 3*254) and (np.amax(r[ignoreColsA-imXRange[0]]) < 3*254):
                ratio5007Hbeta.append(np.amax(rmb[ignoreColsC-imXRange[0]]) / np.amax(rmb[ignoreColsA-imXRange[0]]))
            else:
                ratio5007Hbeta.append(0)
            if (np.amax(r[ignoreColsC-imXRange[0]]) < 3*254) and (np.amax(r[ignoreColsB-imXRange[0]]) < 3*254):
                ratio50074959.append(np.amax(rmb[ignoreColsC-imXRange[0]]) / np.amax(rmb[ignoreColsB-imXRange[0]]))
            else:
                ratio50074959.append(0)
#            print('row = ',row,': 5007 = ',rmb[ignoreColsC-imXRange[0]])
#            print('row = ',row,': Hbeta = ',rmb[ignoreColsA-imXRange[0]])
#            print('row = ',row,': 4959 = ',rmb[ignoreColsB-imXRange[0]])
#            STOP

        plt.plot(np.arange(0,len(ratio5007Hbeta),1)+imYRange[0],ratio5007Hbeta,label='5007/Hbeta')
        #plt.legend()
        plt.xlabel('row')
        plt.ylabel('5007/Hbeta')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_5007-Ratio-Hbeta.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()
        print('len(ratio5007Hbeta) = ',len(ratio5007Hbeta))
        print('mean(5007/Hbeta[',1120-imYRange[0],':',1180-imYRange[0],']) = ',np.mean(ratio5007Hbeta[1120-imYRange[0]:1180-imYRange[0]]))

        plt.plot(np.arange(0,len(ratio50074959),1)+imYRange[0],ratio50074959,label='5007/4959')
        #plt.legend()
        plt.xlabel('row')
        plt.ylabel('5007/4959')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_5007-Ratio-4959.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()
        if iIm == 1:
            print('len(ratio50074959) = ',len(ratio50074959))
            print('mean(5007/4959) for rows 3402-3406 = ',np.mean(ratio50074959[3402-imYRange[0]:3406-imYRange[0]]))
            print('mean(5007/Hbeta) for rows 3402-3406 = ',np.mean(ratio5007Hbeta[3402-imYRange[0]:3406-imYRange[0]]))

        ratio5007HbetaNPix = []
        ratio50074959NPix = []
        print('len(pixelsAboveThreshold5007) = ',len(pixelsAboveThreshold5007))
        for row in np.arange(0,len(pixelsAboveThreshold5007),1):
            print('pixelsAboveThreshold5007[',row,'] = ',len(pixelsAboveThreshold5007[row]),': ',pixelsAboveThreshold5007[row])
            print('pixelsAboveThreshold5007[',row,'][0] = ',len(pixelsAboveThreshold5007[row][0]),': ',pixelsAboveThreshold5007[row][0])
            print('pixelsAboveThreshold4959[',row,'][0] = ',len(pixelsAboveThreshold4959[row][0]),': ',pixelsAboveThreshold4959[row][0])
            print('pixelsAboveThresholdHbeta[',row,'[0]] = ',len(pixelsAboveThresholdHbeta[row][0]),': ',pixelsAboveThresholdHbeta[row][0])
            ratio5007HbetaNPix.append(len(pixelsAboveThreshold5007[row][0]) / len(pixelsAboveThresholdHbeta[row][0]) if len(pixelsAboveThresholdHbeta[row][0])>0 else 0)
            ratio50074959NPix.append(len(pixelsAboveThreshold5007[row][0]) / len(pixelsAboveThreshold4959[row][0]) if len(pixelsAboveThreshold4959[row][0])>0 else 0)

        plt.plot(np.arange(imYRange[0],imYRange[0]+len(pixelsAboveThreshold5007),1),ratio5007HbetaNPix)
        plt.xlabel('row')
        plt.ylabel('nPix(5007)/nPix(Hbeta)')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_nPix5007-Ratio-nPixHbeta.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()

        plt.plot(np.arange(imYRange[0],imYRange[0]+len(pixelsAboveThreshold5007),1),ratio50074959NPix)
        plt.xlabel('row')
        plt.ylabel('nPix(5007)/nPix(4959)')
        plt.title('Plate XII' if iIm == 0 else 'Plate XIII')
        plotName = 'PlateXII' if iIm == 0 else 'PlateXIII'
        plotName += '_nPix5007-Ratio-nPix4959.png'
        plt.savefig(os.path.join(path,plotName))
        plt.show()

def main():
    doFromPdf()

if __name__ == '__main__':
    main()

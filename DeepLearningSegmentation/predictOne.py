import numpy as np
import os
import sys
import model
import cv2


fullname        = os.getcwd()+"/leakyReLu.hdf5"

net = model.unet()
net.load_weights(fullname)


if len(sys.argv) == 3:
    inputFileName= sys.argv[1]
    dir_output=  sys.argv[2]
else :
    print("input not respected")
    print("arg1 : input relief map")
    print("arg2 : output path")

outputFileName=os.path.basename(os.path.splitext(inputFileName)[0])+"SEG.png"

#read image
img = cv2.imread(inputFileName, 0)
h=img.shape[0]
hToPredict=img.shape[0]
w=img.shape[1]
wToPredict=img.shape[1]
#resize to match U-net numfilt
if h%model.numFilt != 0 or w%model.numFilt != 0:
    hToPredict = model.numFilt * round(img.shape[0] / model.numFilt)
    wToPredict = model.numFilt * round(img.shape[1] / model.numFilt)
    img = cv2.resize(img, (wToPredict,hToPredict), interpolation =cv2.INTER_AREA)
#norm
img = img / 255
#tensorflow model take input on 4D
img = np.reshape(img,img.shape+(1,))
img = np.reshape(img,(1,)+img.shape)
#prediction
px 	= net.predict(img, verbose=1)
#extract image
px	= px[0,:,:,0]
#grayscale
px=(px*255).astype(np.uint8)
#resize to original size
if h!=hToPredict or w!=wToPredict :
    px = cv2.resize(px, (w,h), interpolation =cv2.INTER_AREA)
#write
_,px_thresh = cv2.threshold(px,127,255,cv2.THRESH_BINARY)

cv2.imwrite( dir_output+outputFileName, px_thresh );

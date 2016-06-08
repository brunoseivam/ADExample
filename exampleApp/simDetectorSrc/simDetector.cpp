/* simDetector.cpp
 *
 * This is a driver for a simulated area detector.
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  March 20, 2008
 *
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <epicsTime.h>
#include <epicsThread.h>
#include <epicsEvent.h>
#include <epicsMutex.h>
#include <epicsString.h>
#include <epicsStdio.h>
#include <epicsMutex.h>
#include <cantProceed.h>
#include <iocsh.h>

#include <vector>
#include <NDArray.h>

// delete me
#include <pv/ntndarray.h>
using namespace epics::nt;

#include "ADDriver.h"
#include <epicsExport.h>
#include "simDetector.h"

using namespace std;
using namespace epics::pvData;
using namespace epics::pvDatabase;

static const char *driverName = "simDetector";

/* Some systems don't define M_PI in math.h */
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeArray (int sizeX, int sizeY)
{
    switch(SimMode->get()) {
        case SimModeLinearRamp: return computeLinearRampArray<epicsType>(sizeX, sizeY);
        case SimModePeaks:      return computePeaksArray<epicsType>(sizeX, sizeY);
        case SimModeSine:       return computeSineArray<epicsType>(sizeX, sizeY);
    }
    return -1;
}

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeLinearRampArray (int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();
    epicsType *pMono=NULL, *pRed=NULL, *pGreen=NULL, *pBlue=NULL;
    int columnStep=0, rowStep=0;
    epicsType incMono, incRed, incGreen, incBlue;
    int i, j;

    double gain         = ADGain->get();
    double gainX        = SimGainX->get();
    double gainY        = SimGainY->get();
    double gainRed      = SimGainRed->get();
    double gainGreen    = SimGainGreen->get();
    double gainBlue     = SimGainBlue->get();
    int resetImage      = SimResetImage->get();
    int colorMode       = NDColorMode->get();
    double exposureTime = ADAcquireTime->get();

    //gain = 1.0;

    /* The intensity at each pixel[i,j] is:
     * (i * gainX + j* gainY) + imageCounter * gain * exposureTime * 1000. */
    incMono  = (epicsType) (gain      * exposureTime * 1000.);
    incRed   = (epicsType) gainRed   * incMono;
    incGreen = (epicsType) gainGreen * incMono;
    incBlue  = (epicsType) gainBlue  * incMono;

    switch (colorMode) {
    case NDColorModeMono:
        pMono = (epicsType *)rawData;
        break;
    case NDColorModeRGB1:
        columnStep = 3;
        rowStep = 0;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData+1;
        pBlue  = (epicsType *)rawData+2;
        break;
    case NDColorModeRGB2:
        columnStep = 1;
        rowStep = 2 * sizeX;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData + sizeX;
        pBlue  = (epicsType *)rawData + 2*sizeX;
        break;
    case NDColorModeRGB3:
        columnStep = 1;
        rowStep = 0;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData + sizeX*sizeY;
        pBlue  = (epicsType *)rawData + 2*sizeX*sizeY;
        break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);

    if (resetImage) {
        for (i=0; i<sizeY; i++) {
            switch (colorMode) {
            case NDColorModeMono:
                for (j=0; j<sizeX; j++) {
                    (*pMono++) = (epicsType) (incMono * (gainX*j + gainY*i));
                }
                break;
            case NDColorModeRGB1:
            case NDColorModeRGB2:
            case NDColorModeRGB3:
                for (j=0; j<sizeX; j++) {
                    *pRed   = (epicsType) (incRed   * (gainX*j + gainY*i));
                    *pGreen = (epicsType) (incGreen * (gainX*j + gainY*i));
                    *pBlue  = (epicsType) (incBlue  * (gainX*j + gainY*i));
                    pRed   += columnStep;
                    pGreen += columnStep;
                    pBlue  += columnStep;
                }
                pRed   += rowStep;
                pGreen += rowStep;
                pBlue  += rowStep;
                break;
            }
        }
    } else {
        for (i=0; i<sizeY; i++) {
            switch (colorMode) {
            case NDColorModeMono:
                for (j=0; j<sizeX; j++) {
                    *pMono++ += incMono;
                }
                break;
            case NDColorModeRGB1:
            case NDColorModeRGB2:
            case NDColorModeRGB3:
                for (j=0; j<sizeX; j++) {
                    *pRed   += incRed;
                    *pGreen += incGreen;
                    *pBlue  += incBlue;
                    pRed   += columnStep;
                    pGreen += columnStep;
                    pBlue  += columnStep;
                }
                pRed   += rowStep;
                pGreen += rowStep;
                pBlue  += rowStep;
                break;
            }
        }
    }

    pRaw_->setData<epicsType>(data);
    return 0;
}

/** Compute array for array of peaks */
template <typename epicsType> int simDetector::computePeaksArray(int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();
    epicsType *pMono=NULL, *pRed=NULL;
    epicsType *pMono2=NULL, *pRed2=NULL, *pGreen2=NULL, *pBlue2=NULL;
    int columnStep=0;
    int i,j,k,l;
    int minX, maxX, minY,maxY;
    int offsetX, offsetY;
    double gainVariation, noise;
    double gaussX, gaussY;
    double tmpValue;

    int colorMode       = NDColorMode->get();
    double gain         = ADGain->get();
    double gainX        = SimGainX->get();
    double gainY        = SimGainY->get();
    double gainRed      = SimGainRed->get();
    double gainGreen    = SimGainGreen->get();
    double gainBlue     = SimGainBlue->get();
    int peaksStartX     = SimPeakStartX->get();
    int peaksStartY     = SimPeakStartY->get();
    int peaksStepX      = SimPeakStepX->get();
    int peaksStepY      = SimPeakStepY->get();
    int peaksNumX       = SimPeakNumX->get();
    int peaksNumY       = SimPeakNumY->get();
    int peaksWidthX     = SimPeakWidthX->get();
    int peaksWidthY     = SimPeakWidthY->get();
    int peakVariation   = SimPeakHeightVariation->get();
    int noisePct        = SimNoise->get();

    switch (colorMode) {
    case NDColorModeMono:
        pMono = (epicsType *)rawData;
        break;
    case NDColorModeRGB1:
        columnStep = 3;
        pRed   = (epicsType *)rawData;
        break;
    case NDColorModeRGB2:
        columnStep = 1;
        pRed   = (epicsType *)rawData;
        break;
    case NDColorModeRGB3:
        columnStep = 1;
        pRed   = (epicsType *)rawData;
        break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);
    switch (colorMode) {
    case NDColorModeMono:
        // Clear the Image
        pMono2 = pMono;
        for (i = 0; i<sizeY; i++) {
            for (j = 0; j<sizeX; j++) {
                (*pMono2++) = (epicsType)0;
            }
        }
        for (i = 0; i<peaksNumY; i++) {
            for (j = 0; j<peaksNumX; j++) {
                gaussX = 0;
                gaussY = 0;
                if (peakVariation !=0) {
                    gainVariation = 1.0 + (rand()%peakVariation+1)/100.0;
                }
                else{
                    gainVariation = 1.0;
                }
                offsetY = i * peaksStepY + peaksStartY;
                offsetX = j * peaksStepX + peaksStartX;
                minX = (offsetX>4*peaksWidthX) ?(offsetX -4*peaksWidthX):0;
                maxX = (offsetX+4*peaksWidthX<sizeX) ?(offsetX + 4*peaksWidthX):sizeX;
                minY = (offsetY>4*peaksWidthY) ?(offsetY -4*peaksWidthY):0;
                maxY = (offsetY+4*peaksWidthY<sizeY) ?(offsetY + 4*peaksWidthY):sizeY;
                for (k =minY; k<maxY; k++) {
                    pMono2 = pMono + (minX + k*sizeX);
                    for (l=minX; l<maxX; l++) {
                        if (noisePct !=0) {
                            noise = 1.0 + (rand()%noisePct+1)/100.0;
                        }
                        else {
                            noise = 1.0;
                        }
                        gaussY = gainY * exp( -pow((double)(k-offsetY)/(double)peaksWidthY,2.0)/2.0 );
                        gaussX = gainX * exp( -pow((double)(l-offsetX)/(double)peaksWidthX,2.0)/2.0 );
                        tmpValue =  gainVariation*gain * gaussX * gaussY*noise;
                        (*pMono2) += (epicsType)tmpValue;
                        pMono2++;
                    }
                }
            }
        }
        break;
    case NDColorModeRGB1:
    case NDColorModeRGB2:
    case NDColorModeRGB3:
        // Clear the Image
        pRed2 = pRed;
        for (i = 0; i<sizeY; i++) {
            for (j = 0; j<sizeX; j++) {
                (*pRed2++) = (epicsType)0;  //Since we are just clearing the field we will do this with one pointer
                (*pRed2++) = (epicsType)0;
                (*pRed2++) = (epicsType)0;
            }
        }
        for (i = 0; i<peaksNumY; i++) {
            for (j = 0; j<peaksNumX; j++) {
                if (peakVariation !=0) {
                    gainVariation = 1.0 + (rand()%peakVariation+1)/100.0;
                }
                else{
                    gainVariation = 1.0;
                }
                offsetY = i * peaksStepY + peaksStartY;
                offsetX = j * peaksStepX + peaksStartX;
                minX = (offsetX>4*peaksWidthX) ?(offsetX -4*peaksWidthX):0;
                maxX = (offsetX+4*peaksWidthX<sizeX) ?(offsetX + 4*peaksWidthX):sizeX;
                minY = (offsetY>4*peaksWidthY) ?(offsetY -4*peaksWidthY):0;
                maxY = (offsetY+4*peaksWidthY<sizeY) ?(offsetY + 4*peaksWidthY):sizeY;
                for (k =minY; k<maxY; k++) {
                    //Move to the starting point for this peak
                    switch (colorMode) {
                    case NDColorModeRGB1:
                        pRed2 = pRed + (minX*columnStep + k*sizeX*columnStep);
                        pGreen2 = pRed2 + 1;
                        pBlue2 = pRed2 + 2;
                        break;
                    case NDColorModeRGB2:
                        pRed2 = pRed + (minX*columnStep + k*3*sizeX*columnStep);
                        pGreen2 = pRed2 + sizeX;
                        pBlue2 = pRed2 + 2*sizeX;
                        break;
                    case NDColorModeRGB3:
                        pRed2 = pRed + (minX*columnStep + k*sizeX*columnStep);
                        pGreen2 = pRed2 + sizeX*sizeY;
                        pBlue2 = pRed2 + 2*sizeX*sizeY;
                        break;
                    }
                    //Fill in a row for this peak
                    for (l=minX; l<maxX; l++) {
                        if (noisePct !=0) {
                            noise = 1.0 + (rand()%noisePct+1)/100.0;
                        }
                        else {
                            noise = 1.0;
                        }
                        gaussY = gainY * exp( -pow((double)(k-offsetY)/(double)peaksWidthY,2.0)/2.0 );
                        gaussX = gainX * exp( -pow((double)(l-offsetX)/(double)peaksWidthX,2.0)/2.0 );
                        tmpValue =  gainVariation*gain * gaussX * gaussY*noise;
                        (*pRed2) += (epicsType)(gainRed*tmpValue);
                        (*pGreen2) += (epicsType)(gainGreen*tmpValue);
                        (*pBlue2) += (epicsType)(gainBlue*tmpValue);

                        pRed2 += columnStep;
                        pGreen2 += columnStep;
                        pBlue2 += columnStep;
                    }
                }
            }
        }
        break;
    }
    pRaw_->setData<epicsType>(data);
    return 0;
}

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeSineArray(int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();

    epicsType *pMono=NULL, *pRed=NULL, *pGreen=NULL, *pBlue=NULL;
    int columnStep=0, rowStep=0;
    double rndm;
    double xTime, yTime;
    int i, j;

    double gain             = ADGain->get();
    double gainX            = SimGainX->get();
    double gainY            = SimGainY->get();
    double gainRed          = SimGainRed->get();
    double gainGreen        = SimGainGreen->get();
    double gainBlue         = SimGainBlue->get();
    int resetImage          = SimResetImage->get();
    int colorMode           = NDColorMode->get();
    double sineOffset       = SimSineOffset->get();
    double sineNoise        = SimSineNoise->get();
    int xSineOperation      = SimXSineOperation->get();
    double xSine1Amplitude  = SimXSine1Amplitude->get();
    double xSine1Frequency  = SimXSine1Frequency->get();
    double xSine1Phase      = SimXSine1Phase->get();
    double xSine2Amplitude  = SimXSine2Amplitude->get();
    double xSine2Frequency  = SimXSine2Frequency->get();
    double xSine2Phase      = SimXSine2Phase->get();
    int ySineOperation      = SimYSineOperation->get();
    double ySine1Amplitude  = SimYSine1Amplitude->get();
    double ySine1Frequency  = SimYSine1Frequency->get();
    double ySine1Phase      = SimYSine1Phase->get();
    double ySine2Amplitude  = SimYSine2Amplitude->get();
    double ySine2Frequency  = SimYSine2Frequency->get();
    double ySine2Phase      = SimYSine2Phase->get();

    switch (colorMode) {
    case NDColorModeMono:
        pMono = (epicsType *)rawData;
        break;
    case NDColorModeRGB1:
        columnStep = 3;
        rowStep = 0;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData+1;
        pBlue  = (epicsType *)rawData+2;
        break;
    case NDColorModeRGB2:
        columnStep = 1;
        rowStep = 2 * sizeX;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData + sizeX;
        pBlue  = (epicsType *)rawData + 2*sizeX;
        break;
    case NDColorModeRGB3:
        columnStep = 1;
        rowStep = 0;
        pRed   = (epicsType *)rawData;
        pGreen = (epicsType *)rawData + sizeX*sizeY;
        pBlue  = (epicsType *)rawData + 2*sizeX*sizeY;
        break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);

    if (resetImage) {
        if (xSine1_) free(xSine1_);
        if (xSine2_) free(xSine2_);
        if (ySine1_) free(ySine1_);
        if (ySine2_) free(ySine2_);
        xSine1_ = (double *)calloc(sizeX, sizeof(double));
        xSine2_ = (double *)calloc(sizeX, sizeof(double));
        ySine1_ = (double *)calloc(sizeY, sizeof(double));
        ySine2_ = (double *)calloc(sizeY, sizeof(double));
        xSineCounter_ = 0;
        ySineCounter_ = 0;
    }

    for (i=0; i<sizeX; i++) {
        xTime = xSineCounter_++ * gainX / sizeX;
        xSine1_[i] = xSine1Amplitude * sin((xTime  * xSine1Frequency + xSine1Phase/360.) * 2. * M_PI);
        xSine2_[i] = xSine2Amplitude * sin((xTime  * xSine2Frequency + xSine2Phase/360.) * 2. * M_PI);
    }
    for (i=0; i<sizeY; i++) {
        yTime = ySineCounter_++ * gainY / sizeY;
        ySine1_[i] = ySine1Amplitude * sin((yTime  * ySine1Frequency + ySine1Phase/360.) * 2. * M_PI);
        ySine2_[i] = ySine2Amplitude * sin((yTime  * ySine2Frequency + ySine2Phase/360.) * 2. * M_PI);
    }

    if (colorMode == NDColorModeMono) {
        if (xSineOperation == SimSineOperationAdd) {
            for (i=0; i<sizeX; i++) {
                xSine1_[i] = xSine1_[i] + xSine2_[i];
            }
        }
        else {
            for (i=0; i<sizeX; i++) {
                xSine1_[i] = xSine1_[i] * xSine2_[i];
            }
        }
        if (ySineOperation == SimSineOperationAdd) {
            for (i=0; i<sizeY; i++) {
                ySine1_[i] = ySine1_[i] + ySine2_[i];
            }
        }
        else {
            for (i=0; i<sizeY; i++) {
                ySine1_[i] = ySine1_[i] * ySine2_[i];
            }
        }
    }
    for (i=0; i<sizeY; i++) {
        switch (colorMode) {
        case NDColorModeMono:
            for (j=0; j<sizeX; j++) {
                rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                *pMono++ = (epicsType) (gain * (sineOffset + sineNoise*rndm + ySine1_[i] + xSine1_[j]));
            }
            break;
        case NDColorModeRGB1:
        case NDColorModeRGB2:
        case NDColorModeRGB3:
            for (j=0; j<sizeX; j++) {
                rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                *pRed   = (epicsType)(gain * gainRed   * (sineOffset + sineNoise*rndm + xSine1_[j]));
                rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                *pGreen = (epicsType)(gain * gainGreen * (sineOffset + sineNoise*rndm + ySine1_[i]));
                rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                *pBlue  = (epicsType)(gain * gainBlue  * (sineOffset + sineNoise*rndm + (xSine2_[j] + ySine2_[i])/2.));
                pRed   += columnStep;
                pGreen += columnStep;
                pBlue  += columnStep;
            }
            pRed   += rowStep;
            pGreen += rowStep;
            pBlue  += rowStep;
            break;
        }
    }

    pRaw_->setData<epicsType>(data);
    return 0;
}

/** Computes the new image data */
NDArrayPtr simDetector::computeImage()
{
    int xDim=0, yDim=1, colorDim=-1;
    int ndims=0;
    size_t dims[3];
    NDArrayPtr nullImage;
    const char* functionName = "computeImage";
    int status;

    /* NOTE: The caller of this function must have taken the mutex */
    int binX                = ADBinX->get();
    int binY                = ADBinY->get();
    int minX                = ADMinX->get();
    int minY                = ADMinY->get();
    int sizeX               = ADSizeX->get();
    int sizeY               = ADSizeY->get();
    int reverseX            = ADReverseX->get();
    int reverseY            = ADReverseY->get();
    int maxSizeX            = ADMaxSizeX->get();
    int maxSizeY            = ADMaxSizeY->get();
    int colorMode           = NDColorMode->get();
    ScalarType dataType     = (ScalarType)NDDataType->get();
    int resetImage          = SimResetImage->get();

    /* Make sure parameters are consistent, fix them if they are not */
    if (binX < 1) {
        binX = 1;
        ADBinX->put(binX);
    }
    if (binY < 1) {
        binY = 1;
        ADBinY->put(binX);
    }
    if (minX < 0) {
        minX = 0;
        ADMinX->put(minX);
    }
    if (minY < 0) {
        minY = 0;
        ADMinY->put(minY);
    }
    if (minX > maxSizeX-1) {
        minX = maxSizeX-1;
        ADMinX->put(minX);
    }
    if (minY > maxSizeY-1) {
        minY = maxSizeY-1;
        ADMinY->put(minY);
    }
    if (minX+sizeX > maxSizeX) {
        sizeX = maxSizeX-minX;
        ADSizeX->put(sizeX);
    }
    if (minY+sizeY > maxSizeY) {
        sizeY = maxSizeY-minY;
        ADSizeY->put(sizeY);
    }

    switch (colorMode) {
        case NDColorModeMono:
            ndims = 2;
            xDim = 0;
            yDim = 1;
            break;
        case NDColorModeRGB1:
            ndims = 3;
            colorDim = 0;
            xDim     = 1;
            yDim     = 2;
            break;
        case NDColorModeRGB2:
            ndims = 3;
            colorDim = 1;
            xDim     = 0;
            yDim     = 2;
            break;
        case NDColorModeRGB3:
            ndims = 3;
            colorDim = 2;
            xDim     = 0;
            yDim     = 1;
            break;
    }

    if (resetImage) {
        /* Allocate the raw buffer we use to compute images. */
        dims[xDim] = maxSizeX;
        dims[yDim] = maxSizeY;
        if (ndims > 2) dims[colorDim] = 3;

        pRaw_ = mNDArrayPool->alloc(ndims, dims, dataType);

        if (!pRaw_) {
            /*asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
                      "%s:%s: error allocating raw buffer\n",
                      driverName, functionName);*/
            return nullImage;
        }
    }

    switch (dataType) {
        case pvByte:
            status |= computeArray<int8>(maxSizeX, maxSizeY);
            break;
        case pvUByte:
            status |= computeArray<uint8>(maxSizeX, maxSizeY);
            break;
        case pvShort:
            status |= computeArray<int16>(maxSizeX, maxSizeY);
            break;
        case pvUShort:
            status |= computeArray<uint16>(maxSizeX, maxSizeY);
            break;
        case pvInt:
            status |= computeArray<int32>(maxSizeX, maxSizeY);
            break;
        case pvUInt:
            status |= computeArray<uint32>(maxSizeX, maxSizeY);
            break;
        case pvFloat:
            status |= computeArray<float>(maxSizeX, maxSizeY);
            break;
        case pvDouble:
            status |= computeArray<double>(maxSizeX, maxSizeY);
            break;
        default:
            break;
    }

    /* Extract the region of interest with binning.
     * If the entire image is being used (no ROI or binning) that's OK because
     * convertImage detects that case and is very efficient */
    NDDimension_t dimsOut[3];
    NDArray::initDimension(dimsOut[xDim], sizeX, minX, binX, reverseX);
    NDArray::initDimension(dimsOut[yDim], sizeY, minY, binY, reverseY);
    if(ndims > 2) NDArray::initDimension(dimsOut[colorDim], 3);

    vector<NDDimension_t> dimsOutVector;
    for(int i = 0; i < ndims; ++i)
        dimsOutVector.push_back(dimsOut[i]);

    NDArrayPtr result(mNDArrayPool->convert(pRaw_, dataType, dimsOutVector));

    /* We save the most recent image buffer so it can be used in the read() function.
     * Now release it before getting a new version. */
    /*if (this->pArrays[0]) this->pArrays[0]->release();
    status = this->mNDArrayPool->convert(pRaw_,
                                         &this->pArrays[0],
                                         dataType,
                                         dimsOut);
    if (status) {
        //asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
        //            "%s:%s: error allocating buffer in convert()\n",
        //            driverName, functionName);
        return(status);
    }
    //pImage = this->pArrays[0];*/

    //NDArrayInfo_t arrayInfo = pImage->getInfo();

    NDArrayInfo_t arrayInfo = result->getInfo();
    result->setCompressedSize(arrayInfo.totalBytes);
    result->setUncompressedSize(arrayInfo.totalBytes);

    /*NDArraySize->put((int)arrayInfo.totalBytes);
    NDArraySizeX->put((int)pImage->dims[xDim].size);
    NDArraySizeY->put((int)pImage->dims[yDim].size);*/
    SimResetImage->put(0);

    return result;
}

static void simTaskC(void *drvPvt)
{
    simDetector *pPvt = (simDetector *)drvPvt;

    pPvt->simTask();
}

/** This thread calls computeImage to compute new image data and does the callbacks to send it to higher layers.
 * It implements the logic for single, multiple or continuous acquisition. */
void simDetector::simTask()
{
    int status = 0;
    int imageCounter;
    int numImages, numImagesCounter;
    int imageMode;
    int arrayCallbacks;
    int acquire=0;
    NDArrayPtr pImage;
    double acquireTime, acquirePeriod, delay;
    epicsTimeStamp startTime, endTime;
    double elapsedTime;
    //const char *functionName = "simTask";

    this->lock();
    /* Loop forever */
    while (1) {
        /* If we are not acquiring then wait for a semaphore that is given when acquisition is started */
        if (!acquire) {
            /* Release the lock while we wait for an event that says acquire has started, then lock again */
            /*asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                "%s:%s: waiting for acquire to start\n", driverName, functionName);*/
            this->unlock();
            status = epicsEventWait(startEventId_);
            this->lock();
            acquire = 1;
            ADStatusMessage->put("Acquiring data");
            ADNumImagesCounter->put(0);
        }

        /* We are acquiring. */
        /* Get the current time */
        epicsTimeGetCurrent(&startTime);
        imageMode = ADImageMode->get();

        /* Get the exposure parameters */
        acquireTime = ADAcquireTime->get();
        acquirePeriod = ADAcquirePeriod->get();

        ADStatus->put(ADStatusAcquire);

        /* Open the shutter */
        setShutter(ADShutterOpen);

        /* Simulate being busy during the exposure time.  Use epicsEventWaitWithTimeout so that
         * manually stopping the acquisition will work */

        if (acquireTime > 0.0) {
            this->unlock();
            status = epicsEventWaitWithTimeout(stopEventId_, acquireTime);
            this->lock();
        } else {
            status = epicsEventTryWait(stopEventId_);
        }
        if (status == epicsEventWaitOK) {
            acquire = 0;
            if (imageMode == ADImageContinuous) {
                ADStatus->put(ADStatusIdle);
            } else {
                ADStatus->put(ADStatusAborted);
            }
        }

        /* Update the image */
        NDArrayPtr pImage(computeImage());
        if (!pImage) continue;

        /* Close the shutter */
        setShutter(ADShutterClosed);

        if (!acquire) continue;

        ADStatus->put(ADStatusReadout);

        /* Get the current parameters */
        imageCounter = NDArrayCounter->get();
        numImages = ADNumImages->get();
        numImagesCounter = ADNumImagesCounter->get();
        arrayCallbacks = NDArrayCallbacks->get();
        imageCounter++;
        numImagesCounter++;
        NDArrayCounter->put(imageCounter);
        ADNumImagesCounter->put(numImagesCounter);

        /* Put the frame number and time stamp into the buffer */
        pImage->setUniqueId(imageCounter);
        pImage->setTimeStamp(TimeStamp(startTime.secPastEpoch, startTime.nsec));
        //pImage->updateEpicsTimeStamp();

        /* Get any attributes that have been defined for this driver */
        //this->getAttributes(pImage->pAttributeList);

        if (arrayCallbacks) {
            /* Call the NDArray callback */
            /* Must release the lock here, or we can get into a deadlock, because we can
             * block on the plugin lock, and the plugin can be calling us */
            this->unlock();
            /*asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: calling imageData callback\n", driverName, functionName);*/

            NDArrayData->put(pImage->getArray());
            //epicsThreadSleep(1);
            this->lock();
        }

        /* See if acquisition is done */
        if ((imageMode == ADImageSingle) ||
           ((imageMode == ADImageMultiple) &&
            (numImagesCounter >= numImages))) {

            /* First do callback on ADStatus. */
            ADStatusMessage->put("Waiting for acquisition");
            ADStatus->put(ADStatusIdle);

            acquire = 0;
            ADAcquire->put(acquire);
            /*asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: acquisition completed\n", driverName, functionName);*/
        }

        /* If we are acquiring then sleep for the acquire period minus elapsed time. */
        if (acquire) {
            epicsTimeGetCurrent(&endTime);
            elapsedTime = epicsTimeDiffInSeconds(&endTime, &startTime);
            delay = acquirePeriod - elapsedTime;
            /*asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: delay=%f\n",
                    driverName, functionName, delay);*/
            if (delay >= 0.0) {
                /* We set the status to waiting to indicate we are in the period delay */
                ADStatus->put(ADStatusWaiting);
                this->unlock();
                status = epicsEventWaitWithTimeout(stopEventId_, delay);
                this->lock();
                if (status == epicsEventWaitOK) {
                    acquire = 0;
                    if (imageMode == ADImageContinuous) {
                        ADStatus->put(ADStatusIdle);
                    } else {
                        ADStatus->put(ADStatusAborted);
                    }
                }
            }
        }
    }
}
/** Controls the shutter */
void simDetector::setShutter(int open)
{
    if (ADShutterMode->get() == ADShutterModeDetector) {
        /* Simulate a shutter by just changing the status readback */
        ADShutterStatus->put(open);
    } else {
        /* For no shutter or EPICS shutter call the base class method */
        ADDriver::setShutter(open);
    }
}

void simDetector::process (PVRecord const * record)
{
    printf("simDetector process(%s)\n", record->getRecordName().c_str());

    // TODO: recheck if the logic is compatible with original
    if (record == ADAcquire.get()) {
        int value = ADAcquire->get();

        if (value && !acquiring_) {
            ADStatusMessage->put("Acquiring data");
            acquiring_ = true;
            epicsEventSignal(startEventId_);
        }
        if (!value && acquiring_) {
            ADStatusMessage->put("Acquisition stopped");
            if (ADImageMode->get() == ADImageContinuous) {
                ADStatus->put(ADStatusIdle);
            } else {
                ADStatus->put(ADStatusAborted);
            }
            ADStatus->put(ADStatusAcquire);
            epicsEventSignal(stopEventId_);
        }
    } else if ( record == NDDataType.get() ||
                record == NDColorMode.get() ||
                record == SimMode.get() ||
                record == ADAcquireTime.get() ||
                record == ADGain.get() ||
                record == SimGainX.get() ||
                record == SimGainY.get() ||
                record == SimGainRed.get() ||
                record == SimGainGreen.get() ||
                record == SimGainBlue.get()) {
        SimResetImage->put(1);
    } else {
        ADDriver::process(record);
    }

    /*if (status)
        asynPrint(pasynUser, ASYN_TRACE_ERROR,
              "%s:writeInt32 error, status=%d function=%d, value=%d\n",
              driverName, status, function, value);
    else
        asynPrint(pasynUser, ASYN_TRACEIO_DRIVER,
              "%s:writeInt32: function=%d, value=%d\n",
              driverName, function, value);*/
}

/** Report status of the driver.
  * Prints details about the driver if details>0.
  * It then calls the ADDriver::report() method.
  * \param[in] fp File pointed passed by caller where the output is written to.
  * \param[in] details If >0 then driver details are printed.
  */
void simDetector::report(FILE *fp, int details)
{
    fprintf(fp, "Simulation detector\n");
    if (details > 0) {
        fprintf(fp, "  NX, NY:            %d  %d\n", ADSizeX->get(), ADSizeY->get());
        fprintf(fp, "  Data type:         %d\n", NDDataType->get());
    }
    /* Invoke the base class method */
    ADDriver::report(fp, details);
}

/** Constructor for simDetector; most parameters are simply passed to ADDriver::ADDriver.
  * After calling the base class constructor this method creates a thread to compute the simulated detector data,
  * and sets reasonable default values for parameters defined in this class, asynNDArrayDriver and ADDriver.
  * \param[in] portName The name of the asyn port driver to be created.
  * \param[in] maxSizeX The maximum X dimension of the images that this driver can create.
  * \param[in] maxSizeY The maximum Y dimension of the images that this driver can create.
  * \param[in] dataType The initial data type (NDDataType_t) of the images that this driver will create.
  * \param[in] maxBuffers The maximum number of NDArray buffers that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited number of buffers.
  * \param[in] maxMemory The maximum amount of memory that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited amount of memory.
  * \param[in] priority The thread priority for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  * \param[in] stackSize The stack size for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  */
simDetector::simDetector(string const & prefix, int maxSizeX, int maxSizeY,
                         ScalarType dataType, int maxBuffers, size_t maxMemory)

    : ADDriver(prefix, maxBuffers, maxMemory),
      pRaw_(), acquiring_(false), xSine1_(0), xSine2_(0), ySine1_(0), ySine2_(0),
      startEventId_(epicsEventCreate(epicsEventEmpty)),
      stopEventId_(epicsEventCreate(epicsEventEmpty))

{
    //int status = asynSuccess;
    //const char *functionName = "simDetector";

    SimGainX                = createParam<double>(SimGainXString, 1.0);
    SimGainY                = createParam<double>(SimGainYString, 1.0);
    SimGainRed              = createParam<double>(SimGainRedString, 1.0);
    SimGainGreen            = createParam<double>(SimGainGreenString, 1.0);
    SimGainBlue             = createParam<double>(SimGainBlueString, 1.0);
    SimNoise                = createParam<int32>(SimNoiseString, 3);
    SimResetImage           = createParam<int32>(SimResetImageString, 1);
    SimMode                 = createParam<int32>(SimModeString, 0);
    SimPeakStartX           = createParam<int32>(SimPeakStartXString, 1);
    SimPeakStartY           = createParam<int32>(SimPeakStartYString, 1);
    SimPeakWidthX           = createParam<int32>(SimPeakWidthXString, 10);
    SimPeakWidthY           = createParam<int32>(SimPeakWidthYString, 20);
    SimPeakNumX             = createParam<int32>(SimPeakNumXString, 1);
    SimPeakNumY             = createParam<int32>(SimPeakNumYString, 1);
    SimPeakStepX            = createParam<int32>(SimPeakStepXString, 1);
    SimPeakStepY            = createParam<int32>(SimPeakStepYString, 1);
    SimPeakHeightVariation  = createParam<int32>(SimPeakHeightVariationString, 3);
    SimSineOffset           = createParam<double>(SimSineOffsetString);
    SimSineNoise            = createParam<double>(SimSineNoiseString);
    SimXSineOperation       = createParam<int32>(SimXSineOperationString);
    SimXSine1Amplitude      = createParam<double>(SimXSine1AmplitudeString);
    SimXSine1Frequency      = createParam<double>(SimXSine1FrequencyString);
    SimXSine1Phase          = createParam<double>(SimXSine1PhaseString);
    SimXSine2Amplitude      = createParam<double>(SimXSine2AmplitudeString);
    SimXSine2Frequency      = createParam<double>(SimXSine2FrequencyString);
    SimXSine2Phase          = createParam<double>(SimXSine2PhaseString);
    SimYSineOperation       = createParam<int32>(SimYSineOperationString);
    SimYSine1Amplitude      = createParam<double>(SimYSine1AmplitudeString);
    SimYSine1Frequency      = createParam<double>(SimYSine1FrequencyString);
    SimYSine1Phase          = createParam<double>(SimYSine1PhaseString);
    SimYSine2Amplitude      = createParam<double>(SimYSine2AmplitudeString);
    SimYSine2Frequency      = createParam<double>(SimYSine2FrequencyString);
    SimYSine2Phase          = createParam<double>(SimYSine2PhaseString);

    /* Set some default values for parameters */
    ADManufacturer->put("Simulated detector");
    ADModel->put("Basic simulator");
    ADMaxSizeX->put(maxSizeX);
    ADMaxSizeY->put(maxSizeY);
    ADMinX->put(0);
    ADMinY->put(0);
    ADBinX->put(1);
    ADBinY->put(1);
    ADReverseX->put(false);
    ADReverseY->put(false);
    ADSizeX->put(maxSizeX);
    ADSizeY->put(maxSizeY);
    NDArraySizeX->put(maxSizeX);
    NDArraySizeY->put(maxSizeY);
    NDArraySize->put(0);
    NDDataType->put(static_cast<int32>(dataType));
    ADImageMode->put(ADImageContinuous);
    ADAcquireTime->put(.001);
    ADAcquirePeriod->put(.005);
    ADNumImages->put(100);

    int status = (epicsThreadCreate("SimDetTask",
                                epicsThreadPriorityMedium,
                                epicsThreadGetStackSize(epicsThreadStackMedium),
                                (EPICSTHREADFUNC)simTaskC,
                                this) == NULL);
    if (status) {
        printf("simDetector:simDetector epicsThreadCreate failure for image task\n");
        return;
    }
}

/** Configuration command, called directly or from iocsh */
extern "C" int simDetectorConfig(const char *prefix, int maxSizeX, int maxSizeY, int dataType,
                                 int maxBuffers, int maxMemory)
{
    new simDetector(prefix, maxSizeX, maxSizeY, static_cast<ScalarType>(dataType),
                    (maxBuffers < 0) ? 0 : maxBuffers,
                    (maxMemory < 0) ? 0 : maxMemory);
    return 0;
}

/** Code for iocsh registration */
static const iocshArg simDetectorConfigArg0 = {"Prefix", iocshArgString};
static const iocshArg simDetectorConfigArg1 = {"Max X size", iocshArgInt};
static const iocshArg simDetectorConfigArg2 = {"Max Y size", iocshArgInt};
static const iocshArg simDetectorConfigArg3 = {"Data type", iocshArgInt};
static const iocshArg simDetectorConfigArg4 = {"maxBuffers", iocshArgInt};
static const iocshArg simDetectorConfigArg5 = {"maxMemory", iocshArgInt};
static const iocshArg * const simDetectorConfigArgs[] =  {&simDetectorConfigArg0,
                                                          &simDetectorConfigArg1,
                                                          &simDetectorConfigArg2,
                                                          &simDetectorConfigArg3,
                                                          &simDetectorConfigArg4,
                                                          &simDetectorConfigArg5};
static const iocshFuncDef configsimDetector = {"simDetectorConfig", 6, simDetectorConfigArgs};
static void configsimDetectorCallFunc(const iocshArgBuf *args)
{
    simDetectorConfig(args[0].sval, args[1].ival, args[2].ival, args[3].ival,
                      args[4].ival, args[5].ival);
}

static void simDetectorRegister(void)
{
    iocshRegister(&configsimDetector, configsimDetectorCallFunc);
}

extern "C" {
epicsExportRegistrar(simDetectorRegister);
}

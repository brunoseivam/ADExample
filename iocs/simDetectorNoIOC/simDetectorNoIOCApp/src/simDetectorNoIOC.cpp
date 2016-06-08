/* simDetectorNoIOC.cpp
 *
 * This is an example of creating a simDetector and controlling it from outside an IOC
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  October 27, 2014
 *
 */

#include <simDetector.h>
#include <NDPluginDriver.h>
#include <NDPluginROI.h>
#include <pv/pvaClient.h>

#include <pv/pvAccess.h>
#include <pv/channelProviderLocal.h>
#include <pv/serverContext.h>
#include <pv/pvDatabase.h>

using namespace std;
using namespace epics::pvData;
using namespace epics::pvAccess;
using namespace epics::pvaClient;
using namespace epics::pvDatabase;

int main(int argc, char **argv)
{
    cout << "Begin" << endl;

    const string provider("local");
    const string simPrefix("13SIM1:cam1:");
    const string roiPrefix("13SIM1:ROI1:");

    simDetector simD(simPrefix, 1024, 1024, pvByte, 0, 0);
    NDPluginROI roi(roiPrefix, 20, false, "local", simPrefix+string("ArrayData"), 0, 0, 0);

    PvaClientPtr client = PvaClient::get(provider);

    string acquireString(simPrefix + string(ADAcquireString));
    PvaClientChannelPtr acquire = client->channel(acquireString, provider, 1.0);
    PvaClientPutPtr acquirePut = acquire->put();
    PvaClientPutDataPtr acquirePutData = acquirePut->getData();

    string gainString(simPrefix + string(ADGainString));
    PvaClientChannelPtr gain = client->channel(gainString, provider, 1.0);
    PvaClientPutPtr gainPut = gain->put();
    PvaClientPutDataPtr gainPutData = gainPut->getData();

    string callbacksString(simPrefix + string(NDArrayCallbacksString));
    PvaClientChannelPtr arrayCallbacks = client->channel(callbacksString, provider, 1.0);
    PvaClientPutPtr arrayCallbacksPut = arrayCallbacks->put();
    PvaClientPutDataPtr arrayCallbacksPutData = arrayCallbacksPut->getData();

    string roiCallbacksString(roiPrefix + string(NDPluginDriverEnableCallbacksString));
    PvaClientChannelPtr roiEnableCallbacks = client->channel(roiCallbacksString, provider, 1.0);
    PvaClientPutPtr roiEnableCallbacksPut = roiEnableCallbacks->put();
    PvaClientPutDataPtr roiEnableCallbacksPutData = roiEnableCallbacksPut->getData();

    string roiDataTypeString(roiPrefix + string(NDPluginROIDataTypeString));
    PvaClientChannelPtr roiDataType = client->channel(roiDataTypeString, provider, 1.0);
    PvaClientPutPtr roiDataTypePut = roiDataType->put();
    PvaClientPutDataPtr roiDataTypePutData = roiDataTypePut->getData();

    string roiDim0MinString(roiPrefix + string(NDPluginROIDim0MinString));
    PvaClientChannelPtr roiDim0Min = client->channel(roiDim0MinString, provider, 1.0);
    PvaClientPutPtr roiDim0MinPut = roiDim0Min->put();
    PvaClientPutDataPtr roiDim0MinPutData = roiDim0MinPut->getData();

    string roiBlockingCallbacksString(roiPrefix + string(NDPluginDriverBlockingCallbacksString));
    PvaClientChannelPtr roiBlockingCallbacks = client->channel(roiBlockingCallbacksString, provider, 1.0);
    PvaClientPutPtr roiBlockingCallbacksPut = roiBlockingCallbacks->put();
    PvaClientPutDataPtr roiBlockingCallbacksPutData = roiBlockingCallbacksPut->getData();

    //roiDim0MinPutData->putDouble(499);
    //roiDim0MinPut->put();

    roiDataTypePutData->putDouble(1.0);
    roiDataTypePut->put();

    roiEnableCallbacksPutData->putDouble(1.0);
    roiEnableCallbacksPut->put();

    roiBlockingCallbacksPutData->putString("true");
    roiBlockingCallbacksPut->put();

    gainPutData->putDouble(1.0);
    gainPut->put();

    arrayCallbacksPutData->putDouble(1.0);
    arrayCallbacksPut->put();

    acquirePutData->putDouble(1);
    acquirePut->put();

    startPVAServer("local", 0, false, true);

    for(;;);
    cout << "End" << endl;

  return 0;
}

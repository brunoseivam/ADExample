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

#include <pv/pvPortClient.h>
#include <pv/scalarClient.h>

using namespace std;
using namespace epics::pvData;
using namespace epics::pvAccess;
using namespace epics::pvPortClient;

int main(int argc, char **argv)
{
    const string provider("local");
    const string simPrefix("13SIM1:cam1:");
    const string roiPrefix("13SIM1:ROI1:");

    simDetector simD(simPrefix, 1024, 1024, pvByte, 0, 0);
    NDPluginROI roi(roiPrefix, 20, false, provider, simPrefix+string(NDArrayDataString), 0, 0, 0);

    IntParamClientPtr acquire(IntParamClient::create(provider, simPrefix + string(ADAcquireString)));
    DoubleParamClientPtr gain(DoubleParamClient::create(provider, simPrefix + string(ADGainString)));
    IntParamClientPtr arrayCallbacks(IntParamClient::create(provider, simPrefix + string(NDArrayCallbacksString)));

    IntParamClientPtr roiEnableCallbacks(IntParamClient::create(provider, roiPrefix + string(NDPluginDriverEnableCallbacksString)));
    IntParamClientPtr roiDataType(IntParamClient::create(provider, roiPrefix + string(NDPluginROIDataTypeString)));
    BooleanParamClientPtr roiBlockingCallbacks(BooleanParamClient::create(provider, roiPrefix + string(NDPluginDriverBlockingCallbacksString)));

    roiDataType->put(1);
    roiEnableCallbacks->put(1);
    roiBlockingCallbacks->put(true);

    gain->put(1.0);
    arrayCallbacks->put(1);
    acquire->put(1);

    startPVAServer("local", 0, false, true);

    return 0;
}

#include <epicsEvent.h>
#include <ADDriver.h>
#include <NDArray.h>

/** Simulation detector driver; demonstrates most of the features that areaDetector drivers can support. */
class epicsShareClass simDetector : public ADDriver {
public:
    simDetector(std::string const & prefix, int maxSizeX, int maxSizeY,
                epics::pvData::ScalarType dataType, int maxBuffers, size_t maxMemory);

    /* These are the methods that we override from ADDriver */
    void process (epics::pvDatabase::PVRecord const * record);
    virtual void setShutter(int open);
    virtual void report(FILE *fp, int details);

    void simTask(); /**< Should be private, but gets called from C, so must be public */

protected:
    epics::pvPortDriver::DoubleParamPtr SimGainX;
    epics::pvPortDriver::DoubleParamPtr SimGainY;
    epics::pvPortDriver::DoubleParamPtr SimGainRed;
    epics::pvPortDriver::DoubleParamPtr SimGainGreen;
    epics::pvPortDriver::DoubleParamPtr SimGainBlue;
    epics::pvPortDriver::IntParamPtr    SimNoise;
    epics::pvPortDriver::IntParamPtr    SimResetImage;
    epics::pvPortDriver::IntParamPtr    SimMode;
    epics::pvPortDriver::IntParamPtr    SimPeakStartX;
    epics::pvPortDriver::IntParamPtr    SimPeakStartY;
    epics::pvPortDriver::IntParamPtr    SimPeakWidthX;
    epics::pvPortDriver::IntParamPtr    SimPeakWidthY;
    epics::pvPortDriver::IntParamPtr    SimPeakNumX;
    epics::pvPortDriver::IntParamPtr    SimPeakNumY;
    epics::pvPortDriver::IntParamPtr    SimPeakStepX;
    epics::pvPortDriver::IntParamPtr    SimPeakStepY;
    epics::pvPortDriver::IntParamPtr    SimPeakHeightVariation;
    epics::pvPortDriver::DoubleParamPtr SimSineOffset;
    epics::pvPortDriver::DoubleParamPtr SimSineNoise;
    epics::pvPortDriver::IntParamPtr    SimXSineOperation;
    epics::pvPortDriver::DoubleParamPtr SimXSine1Amplitude;
    epics::pvPortDriver::DoubleParamPtr SimXSine1Frequency;
    epics::pvPortDriver::DoubleParamPtr SimXSine1Phase;
    epics::pvPortDriver::DoubleParamPtr SimXSine2Amplitude;
    epics::pvPortDriver::DoubleParamPtr SimXSine2Frequency;
    epics::pvPortDriver::DoubleParamPtr SimXSine2Phase;
    epics::pvPortDriver::IntParamPtr    SimYSineOperation;
    epics::pvPortDriver::DoubleParamPtr SimYSine1Amplitude;
    epics::pvPortDriver::DoubleParamPtr SimYSine1Frequency;
    epics::pvPortDriver::DoubleParamPtr SimYSine1Phase;
    epics::pvPortDriver::DoubleParamPtr SimYSine2Amplitude;
    epics::pvPortDriver::DoubleParamPtr SimYSine2Frequency;
    epics::pvPortDriver::DoubleParamPtr SimYSine2Phase;

private:
    /* These are the methods that are new to this class */
    template <typename epicsType> int computeArray(int sizeX, int sizeY);
    template <typename epicsType> int computeLinearRampArray(int sizeX, int sizeY);
    template <typename epicsType> int computePeaksArray(int sizeX, int sizeY);
    template <typename epicsType> int computeSineArray(int sizeX, int sizeY);
    NDArrayPtr computeImage();

    /* Our data */
    epicsEventId startEventId_;
    epicsEventId stopEventId_;
    NDArrayPtr pRaw_;
    bool acquiring_;
    double *xSine1_;
    double *xSine2_;
    double *ySine1_;
    double *ySine2_;
    double xSineCounter_;
    double ySineCounter_;
};

typedef enum {
    SimModeLinearRamp,
    SimModePeaks,
    SimModeSine
} SimModes_t;

typedef enum {
    SimSineOperationAdd,
    SimSineOperationMultiply
} SimSineOperation_t;

#define SimGainXString                "SIM_GAIN_X"
#define SimGainYString                "SIM_GAIN_Y"
#define SimGainRedString              "SIM_GAIN_RED"
#define SimGainGreenString            "SIM_GAIN_GREEN"
#define SimGainBlueString             "SIM_GAIN_BLUE"
#define SimNoiseString                "SIM_NOISE"
#define SimResetImageString           "RESET_IMAGE"
#define SimModeString                 "SIM_MODE"
#define SimPeakStartXString           "SIM_PEAK_START_X"
#define SimPeakStartYString           "SIM_PEAK_START_Y"
#define SimPeakWidthXString           "SIM_PEAK_WIDTH_X"
#define SimPeakWidthYString           "SIM_PEAK_WIDTH_Y"
#define SimPeakNumXString             "SIM_PEAK_NUM_X"
#define SimPeakNumYString             "SIM_PEAK_NUM_Y"
#define SimPeakStepXString            "SIM_PEAK_STEP_X"
#define SimPeakStepYString            "SIM_PEAK_STEP_Y"
#define SimPeakHeightVariationString  "SIM_PEAK_HEIGHT_VARIATION"
#define SimSineOffsetString           "SIM_SINE_OFFSET"
#define SimSineNoiseString            "SIM_SINE_NOISE"
#define SimXSineOperationString       "SIM_XSINE_OPERATION"
#define SimXSine1AmplitudeString      "SIM_XSINE1_AMPLITUDE"
#define SimXSine1FrequencyString      "SIM_XSINE1_FREQUENCY"
#define SimXSine1PhaseString          "SIM_XSINE1_PHASE"
#define SimXSine2AmplitudeString      "SIM_XSINE2_AMPLITUDE"
#define SimXSine2FrequencyString      "SIM_XSINE2_FREQUENCY"
#define SimXSine2PhaseString          "SIM_XSINE2_PHASE"
#define SimYSineOperationString       "SIM_YSINE_OPERATION"
#define SimYSine1AmplitudeString      "SIM_YSINE1_AMPLITUDE"
#define SimYSine1FrequencyString      "SIM_YSINE1_FREQUENCY"
#define SimYSine1PhaseString          "SIM_YSINE1_PHASE"
#define SimYSine2AmplitudeString      "SIM_YSINE2_AMPLITUDE"
#define SimYSine2FrequencyString      "SIM_YSINE2_FREQUENCY"
#define SimYSine2PhaseString          "SIM_YSINE2_PHASE"


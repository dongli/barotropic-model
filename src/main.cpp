#include "BarotropicModel.h"
#include "RossbyHaurwitzTestCase.h"

int main(int argc, const char *argv[])
{
    BarotropicModel model;
    RossbyHaurwitzTestCase testCase;

    TimeManager timeManager;
    Time startTime, endTime(86400);

    timeManager.init(startTime, endTime, 360);

    model.init(80, 41);
    testCase.calcInitCond(model);

    model.run(timeManager);

    return 0;
}

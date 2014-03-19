#include "BarotropicModel.h"
#include "RossbyHaurwitzTestCase.h"

int main(int argc, const char *argv[])
{
    BarotropicModel model;
    RossbyHaurwitzTestCase testCase;

    TimeManager timeManager;
    Time startTime, endTime(68*DAYS);

    timeManager.init(startTime, endTime, 4*MINUTES);

    model.init(80, 41);
    testCase.calcInitCond(model);

    model.run(timeManager);

    return 0;
}

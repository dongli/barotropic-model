#include "barotropic_model.h"

using namespace barotropic_model;

int main(int argc, const char *argv[])
{
    BarotropicModel_A_ImplicitMidpoint model;
    RossbyHaurwitzTestCase testCase;

    TimeManager timeManager;
    ptime startTime(date(2000, 1, 1));
    ptime endTime = startTime+days(120);

    timeManager.init(startTime, endTime, minutes(4));

    model.init(timeManager, 80, 41);
    testCase.calcInitCond(model);

    model.run();

    return 0;
}

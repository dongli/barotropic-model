#include "BarotropicModel_C_ImplicitMidpoint.h"

namespace barotropic_model {

BarotropicModel_C_ImplicitMidpoint::BarotropicModel_C_ImplicitMidpoint() {
    REPORT_ONLINE;
}

BarotropicModel_C_ImplicitMidpoint::~BarotropicModel_C_ImplicitMidpoint() {
    REPORT_OFFLINE;
}

void BarotropicModel_C_ImplicitMidpoint::init(TimeManager &timeManager,
                                              int numLon, int numLat) {
    this->timeManager = &timeManager;
    // initialize IO manager
    io.init(timeManager);
    // initialize domain
    domain = new Domain(2);
    domain->setRadius(6.371e6);
    // initialize mesh
    mesh = new Mesh(*domain);
    mesh->init(numLon, numLat);
    // create variables
    u.create("u", "m s-1", "zonal wind speed", *mesh, X_FACE, HAS_HALF_LEVEL);
    v.create("v", "m s-1", "meridional wind speed", *mesh, Y_FACE, HAS_HALF_LEVEL);
    gd.create("gd", "m", "geopotential height", *mesh, CENTER, HAS_HALF_LEVEL);
    ut.create("ut", "(m s-1)2", "transformed zonal wind speed", *mesh, X_FACE, HAS_HALF_LEVEL);
    vt.create("vt", "(m s-1)2", "transformed meridional wind speed", *mesh, Y_FACE, HAS_HALF_LEVEL);
    gdt.create("gdt", "m2 s-2", "transformed geopotential height", *mesh, CENTER, HAS_HALF_LEVEL);
    dut.create("dut", "m s-2", "zonal wind speed tendency", *mesh, X_FACE);
    dvt.create("dvt", "m s-2", "meridional zonal speed tendency", *mesh, Y_FACE);
    dgd.create("dgd", "m-2 s-1", "geopotential height tendency", *mesh, CENTER);
    gdu.create("gdu", "m2 s-1", "ut * gdt", *mesh, X_FACE);
    gdv.create("gdv", "m2 s-1", "vt * gdt", *mesh, Y_FACE);
    // set coefficients
}

void BarotropicModel_C_ImplicitMidpoint::run() {
    int fileIdx = io.registerOutputFile(*mesh, "output", IOFrequencyUnit::DAYS, 0.5);
    io.file(fileIdx).registerOutputField<double, 2, FULL_DIMENSION>(3, &u, &v, &gd);
    // -------------------------------------------------------------------------
    // output initial condition
    io.create(fileIdx);
    io.output<double>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
    io.close(fileIdx);
    // -------------------------------------------------------------------------
    // main integration loop
    while (!timeManager->isFinished()) {
        integrate(oldTimeIdx, timeManager->getStepSize());
        timeManager->advance();
        oldTimeIdx.shift();
        io.create(fileIdx);
        io.output<double>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
        io.close(fileIdx);
    }
}

void BarotropicModel_C_ImplicitMidpoint::integrate(const TimeLevelIndex &oldTimeIdx,
                                                   double dt) {
    // -------------------------------------------------------------------------
    // set time level indices
    halfTimeIdx = oldTimeIdx+0.5;
    newTimeIdx = oldTimeIdx+1;
}

double BarotropicModel_C_ImplicitMidpoint::calcTotalEnergy(const TimeLevelIndex &timeIdx) const {
    double totalEnergy = 0;
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, HALF); ++i) {
            totalEnergy += pow(ut(timeIdx, i, j), 2)*cosLatFull[j];
        }
    }
    for (int j = 0; j < mesh->getNumGrid(1, HALF); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            totalEnergy += pow(vt(timeIdx, i, j), 2)*cosLatHalf[j];
        }
    }
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            totalEnergy += pow(gd(timeIdx, i, j)+ghs(i, j), 2)*cosLatFull[j];
        }
    }
    return totalEnergy;
}

double BarotropicModel_C_ImplicitMidpoint::calcTotalMass(const TimeLevelIndex &timeIdx) const {
    double totalMass = 0;
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            totalMass += gd(timeIdx, i, j)*cosLatFull[j];
        }
    }
    return totalMass;
}

}

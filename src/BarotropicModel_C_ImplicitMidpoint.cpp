#include "BarotropicModel_C_ImplicitMidpoint.h"

namespace barotropic_model {

BarotropicModel_C_ImplicitMidpoint::BarotropicModel_C_ImplicitMidpoint() {
    REPORT_ONLINE;
}

BarotropicModel_C_ImplicitMidpoint::~BarotropicModel_C_ImplicitMidpoint() {
    REPORT_OFFLINE;
}

void BarotropicModel_C_ImplicitMidpoint::init(int numLon, int numLat) {
    // -------------------------------------------------------------------------
    // initialize domain
    domain = new Domain(2);
    domain->setRadius(6.371e6);
    // -------------------------------------------------------------------------
    // initialize mesh
    mesh = new Mesh(*domain);
    mesh->init(numLon, numLat);
    // -------------------------------------------------------------------------
    // create variables
    u.create("u", "m s-1", "zonal wind speed", *mesh, X_FACE, HAS_HALF_LEVEL);
    v.create("v", "m s-1", "meridional wind speed", *mesh, Y_FACE, HAS_HALF_LEVEL);
    gd.create("gd", "m", "geopotential height", *mesh, CENTER, HAS_HALF_LEVEL);

    ut.create("ut", "(m s-1)2", "transformed zonal wind speed", *mesh, X_FACE, HAS_HALF_LEVEL);
    vt.create("vt", "(m s-1)2", "transformed meridional wind speed", *mesh, Y_FACE, HAS_HALF_LEVEL);
    ght.create("ght", "m2 s-2", "transformed geopotential height", *mesh, CENTER, HAS_HALF_LEVEL);

    dut.create("dut", "m s-2", "zonal wind speed tendency", *mesh, X_FACE);
    dvt.create("dvt", "m s-2", "meridional zonal speed tendency", *mesh, Y_FACE);
    dgd.create("dgd", "m-2 s-1", "geopotential height tendency", *mesh, CENTER);

    ghu.create("ghu", "m2 s-1", "ut * ght", *mesh, X_FACE);
    ghv.create("gdv", "m2 s-1", "vt * ght", *mesh, Y_FACE);

//    uut.create("uut", "(m s-1)3", "u * ut", *mesh, X_FACE);
//    vut.create("vut", "(m s-1)3", "v * ut", *mesh, XY_VERTEX);
//    
//    uvt.create("uvt", "(m s-1)3", "u * vt", *mesh, XY_VERTEX);
//    vvt.create("vvt", "(m s-1)3", "v * vt", *mesh, Y_FACE);
    // -------------------------------------------------------------------------
    // set coefficients
}

void BarotropicModel_C_ImplicitMidpoint::run(TimeManager &timeManager) {
    // -------------------------------------------------------------------------
    // initialize IO manager
    io.init(timeManager);
    int fileIdx = io.registerOutputFile(*mesh, "output", IOFrequencyUnit::DAYS, 0.5);
    io.file(fileIdx).registerOutputField<double, 2, FULL_DIMENSION>(3, &u, &v, &gd);
    // -------------------------------------------------------------------------
    // output initial condition
    io.create(fileIdx);
    io.output<double>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
    io.close(fileIdx);
    // -------------------------------------------------------------------------
    // main integration loop
    while (!timeManager.isFinished()) {
        integrate(oldTimeIdx, timeManager.getStepSize());
        timeManager.advance();
        oldTimeIdx.shift();
        io.create(fileIdx);
        io.output<double>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
        io.close(fileIdx);
    }
}

void BarotropicModel_C_ImplicitMidpoint::integrate(const TimeLevelIndex &oldTimeIdx,
                                                   double dt) {
    
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
            totalEnergy += pow(gd(timeIdx, i, j), 2)*cosLatFull[j];
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

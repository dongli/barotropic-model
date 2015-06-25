#include "BarotropicModel_A_ImplicitMidpoint.h"

namespace barotropic_model {

BarotropicModel_A_ImplicitMidpoint::BarotropicModel_A_ImplicitMidpoint() {
    REPORT_ONLINE;
}

BarotropicModel_A_ImplicitMidpoint::~BarotropicModel_A_ImplicitMidpoint() {
    REPORT_OFFLINE;
}

void BarotropicModel_A_ImplicitMidpoint::
init(TimeManager &timeManager, int numLon, int numLat) {
    this->timeManager = &timeManager;
    // Initialize the IO manager.
    io.init(timeManager);
    // Initialize the domain.
    _domain = new Domain(2);
    domain().radius() = 6.371e6;
    // Initialize the mesh.
    _mesh = new Mesh(domain());
    mesh().init(numLon, numLat);
    dlon = mesh().gridInterval(0, FULL, 0);
    dlat = mesh().gridInterval(1, FULL, 0); // Assume the equidistance grids.
    // Create the variables.
    u.create("u", "m s-1", "zonal wind speed", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    v.create("v", "m s-1", "meridional wind speed", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    gd.create("gd", "m2 s-2", "geopotential depth", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    ghs.create("ghs", "m2 s-2", "surface geopotential", mesh(), CENTER, 2);
    ut.create("ut", "(m s-1)*m-2", "transformed zonal wind speed", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    vt.create("vt", "(m s-1)*m-2", "transformed meridional wind speed", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    gdt.create("gdt", "m-2", "transformed geopotential height", mesh(), CENTER, 2, HAS_HALF_LEVEL);
    dut.create("dut", "m s-2", "zonal wind speed tendency", mesh(), CENTER, 2);
    dvt.create("dvt", "m s-2", "meridional zonal speed tendency", mesh(), CENTER, 2);
    dgd.create("dgd", "m-2 s-1", "geopotential depth tendency", mesh(), CENTER, 2);
    gdu.create("gdu", "m2 s-1", "ut * gdt", mesh(), CENTER, 2);
    gdv.create("gdv", "m2 s-1", "vt * gdt", mesh(), CENTER, 2);
    fu.create("fu", "* m s-1", "* * u", mesh(), CENTER, 2);
    fv.create("fv", "* m s-1", "* * v", mesh(), CENTER, 2);
    // Set some coefficients.
    // Note: Some coefficients containing cos(lat) will be specialized at Poles.
    cosLat.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        cosLat[j] = mesh().cosLat(FULL, j);
    }
    cosLat[mesh().js(FULL)] = mesh().cosLat(HALF, mesh().js(HALF))*0.25;
    cosLat[mesh().je(FULL)] = mesh().cosLat(HALF, mesh().je(HALF))*0.25;
    tanLat.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        tanLat[j] = mesh().tanLat(FULL, j);
    }
    tanLat[mesh().js(FULL)] = -1/cosLat[mesh().js(FULL)];
    tanLat[mesh().je(FULL)] =  1/cosLat[mesh().je(FULL)];
    factorCor.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        factorCor[j] = 2*OMEGA*mesh().sinLat(FULL, j);
    }
    factorCur.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        factorCur[j] = tanLat[j]/domain().radius();
    }
    factorLon.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        factorLon[j] = 1/(2*dlon*domain().radius()*cosLat[j]);
    }
    factorLat.set_size(mesh().numGrid(1, FULL));
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        factorLat[j] = 1/(2*dlat*domain().radius()*cosLat[j]);
    }
    // Set the variables on the Poles.
    for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
        dut(i, mesh().js(FULL)) = 0.0; dut(i, mesh().je(FULL)) = 0.0;
        dvt(i, mesh().js(FULL)) = 0.0; dvt(i, mesh().je(FULL)) = 0.0;
        gdu(i, mesh().js(FULL)) = 0.0; gdu(i, mesh().je(FULL)) = 0.0;
        gdv(i, mesh().js(FULL)) = 0.0; gdv(i, mesh().je(FULL)) = 0.0;
    }
} // init

void BarotropicModel_A_ImplicitMidpoint::
input(const string &fileName) {
    int fileIdx = io.addInputFile(mesh(), fileName);
    io.file(fileIdx).addField("double", FULL_DIMENSION, {&u, &v, &gd});
    io.file(fileIdx).addField("double", FULL_DIMENSION, {&ghs});
    io.open(fileIdx);
    io.updateTime(fileIdx, *timeManager);
    io.input<double, 2>(fileIdx, oldTimeIdx, {&u, &v, &gd});
    io.input<double>(fileIdx, {&ghs});
    io.close(fileIdx);
    io.removeFile(fileIdx);
    u.applyBndCond(oldTimeIdx);
    v.applyBndCond(oldTimeIdx);
    gd.applyBndCond(oldTimeIdx);
    ghs.applyBndCond();
} // input

void BarotropicModel_A_ImplicitMidpoint::
run() {
    // Add the output fields.
    StampString filePattern("output.%5s.nc");
    int fileIdx = io.addOutputFile(mesh(), filePattern, hours(1));
    io.file(fileIdx).addField("double", FULL_DIMENSION, {&u, &v, &gd});
    io.file(fileIdx).addField("double", FULL_DIMENSION, {&ghs});
    // Output the initial condition.
    io.create(fileIdx);
    io.output<double, 2>(fileIdx, oldTimeIdx, {&u, &v, &gd});
    io.output<double>(fileIdx, {&ghs});
    io.close(fileIdx);
    // Start the main integration loop.
    while (!timeManager->isFinished()) {
        integrate(oldTimeIdx, timeManager->stepSizeInSeconds());
        timeManager->advance();
        oldTimeIdx.shift();
        io.create(fileIdx);
        io.output<double, 2>(fileIdx, oldTimeIdx, {&u, &v, &gd});
        io.output<double>(fileIdx, {&ghs});
        io.close(fileIdx);
    }
} // run

void BarotropicModel_A_ImplicitMidpoint::
integrate(const TimeLevelIndex<2> &oldTimeIdx, double dt) {
    // Set time level indices.
    halfTimeIdx = oldTimeIdx+0.5;
    newTimeIdx = oldTimeIdx+1;
    // Copy the old variables to the new ones to start the iteration.
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
            u(newTimeIdx, i, j) = u(oldTimeIdx, i, j);
            v(newTimeIdx, i, j) = v(oldTimeIdx, i, j);
            gd(newTimeIdx, i, j) = gd(oldTimeIdx, i, j);
            gdt(oldTimeIdx, i, j) = sqrt(gd(oldTimeIdx, i, j));
            gdt(newTimeIdx, i, j) = gdt(oldTimeIdx, i, j);
            ut(oldTimeIdx, i, j) = u(oldTimeIdx, i, j)*gdt(oldTimeIdx, i, j);
            ut(newTimeIdx, i, j) = ut(oldTimeIdx, i, j);
            vt(oldTimeIdx, i, j) = v(oldTimeIdx, i, j)*gdt(oldTimeIdx, i, j);
            vt(newTimeIdx, i, j) = vt(oldTimeIdx, i, j);
        }
    }
    // Get the old total energy and mass.
    double e0 = calcTotalEnergy(oldTimeIdx);
    double m0 = calcTotalMass(oldTimeIdx);
#ifndef NDEBUG
    cout << "iteration ";
#endif
    cout << "energy: ";
    cout << std::fixed << setw(20) << setprecision(2) << e0 << "  ";
    cout << "mass: ";
    cout << setw(20) << setprecision(2) << m0 << endl;
    // Run iterations.
    for (int iter = 1; iter <= 8; ++iter) {
        // Calculate the variables on the half time step.
        for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
            for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
                u(halfTimeIdx, i, j) = (u(oldTimeIdx, i, j)+u(newTimeIdx, i, j))*0.5;
                v(halfTimeIdx, i, j) = (v(oldTimeIdx, i, j)+v(newTimeIdx, i, j))*0.5;
                gd(halfTimeIdx, i, j) = (gd(oldTimeIdx, i, j)+gd(newTimeIdx, i, j))*0.5;
                ut(halfTimeIdx, i, j) = (ut(oldTimeIdx, i, j)+ut(newTimeIdx, i, j))*0.5;
                vt(halfTimeIdx, i, j) = (vt(oldTimeIdx, i, j)+vt(newTimeIdx, i, j))*0.5;
                gdt(halfTimeIdx, i, j) = (gdt(oldTimeIdx, i, j)+gdt(newTimeIdx, i, j))*0.5;
            }
        }
        // Update the geopotential height.
        calcGeopotentialDepthTendency(halfTimeIdx);
        for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
            for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
                gd(newTimeIdx, i, j) = gd(oldTimeIdx, i, j)-dt*dgd(i, j);
            }
        }
        gd.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        // Transform the geopotential height.
        for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
            for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
                gdt(newTimeIdx, i, j) = sqrt(gd(newTimeIdx, i, j));
            }
        }
        gdt.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        // Update the velocity.
        calcZonalWindTendency(halfTimeIdx);
        calcMeridionalWindTendency(halfTimeIdx);
        for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
            for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
                ut(newTimeIdx, i, j) = ut(oldTimeIdx, i, j)-dt*dut(i, j);
                vt(newTimeIdx, i, j) = vt(oldTimeIdx, i, j)-dt*dvt(i, j);
            }
        }
        ut.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        vt.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        // Transform back velocity on the new time level.
        for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
            for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
                u(newTimeIdx, i, j) = ut(newTimeIdx, i, j)/gdt(newTimeIdx, i, j);
                v(newTimeIdx, i, j) = vt(newTimeIdx, i, j)/gdt(newTimeIdx, i, j);
            }
        }
        u.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        v.applyBndCond(newTimeIdx, UPDATE_HALF_LEVEL);
        // Get the new total energy and mass.
        double e1 = calcTotalEnergy(newTimeIdx);
        // TODO: Figure out how this early iteration abortion works.
        if (fabs(e1-e0)*2/(e1+e0) < 5.0e-15) {
            break;
        }
#ifndef NDEBUG
        double m1 = calcTotalMass(newTimeIdx);
        cout << setw(9) << iter;
        cout << " energy: ";
        cout << std::fixed << setw(20) << setprecision(2) << e1 << "  ";
        cout << "mass: ";
        cout << setw(20) << setprecision(2) << m1 << " ";
        cout << "energy bias: ";
        cout << setw(20) << setprecision(16) << fabs(e1-e0)*2/(e1+e0) << endl;
#endif
    }
} // integrate

double BarotropicModel_A_ImplicitMidpoint::
calcTotalEnergy(const TimeLevelIndex<2> &timeIdx) const {
    double totalEnergy = 0.0;
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            totalEnergy += (pow(ut(timeIdx, i, j), 2)+
                            pow(vt(timeIdx, i, j), 2)+
                            pow(gd(timeIdx, i, j)+ghs(i, j), 2))*cosLat[j];
        }
    }
    return totalEnergy;
} // calcTotalEnergy

double BarotropicModel_A_ImplicitMidpoint::
calcTotalMass(const TimeLevelIndex<2> &timeIdx) const {
    double totalMass = 0.0;
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            totalMass += gd(timeIdx, i, j)*cosLat[j];
        }
    }
    return totalMass;
} // calcTotalMass

/**
 *  Input: ut, vt, gdt
 *  Intermediate: gdu, gdv
 *  Output: dgd
 */
void BarotropicModel_A_ImplicitMidpoint::
calcGeopotentialDepthTendency(const TimeLevelIndex<2> &timeIdx) {
    // calculate intermediate variables
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
            gdu(i, j) = ut(timeIdx, i, j)*gdt(timeIdx, i, j);
            gdv(i, j) = vt(timeIdx, i, j)*gdt(timeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            dgd(i, j) = (gdu(i+1, j)-gdu(i-1, j))*factorLon[j]+
                        (gdv(i, j+1)-gdv(i, j-1))*factorLat[j];
        }
    }
    // pole grids
    // last character 's' and 'n' mean 'Sorth Pole' and 'North Pole' respectively
    int js = mesh().js(FULL), jn = mesh().je(FULL);
    double dgds = 0.0, dgdn = 0.0;
    for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
        dgds += gdv(i, js+1);
        dgdn -= gdv(i, jn-1);
    }
    dgds *= factorLat[js]/mesh().numGrid(0, FULL);
    dgdn *= factorLat[jn]/mesh().numGrid(0, FULL);
    for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
        dgd(i, js) = dgds;
        dgd(i, jn) = dgdn;
    }
#ifndef NDEBUG
    double tmp = 0.0;
    for (int j = mesh().js(FULL); j <= mesh().je(FULL); ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            tmp += dgd(i, j)*cosLat[j];
        }
    }
    assert(fabs(tmp) < 1.0e-10);
#endif
} // calcGeopotentialDepthTendency

void BarotropicModel_A_ImplicitMidpoint::
calcZonalWindTendency(const TimeLevelIndex<2> &timeIdx) {
    calcZonalWindAdvection(timeIdx);
    calcZonalWindCoriolis(timeIdx);
    calcZonalWindPressureGradient(timeIdx);
} // calcZonalWindTendency

void BarotropicModel_A_ImplicitMidpoint::
calcMeridionalWindTendency(const TimeLevelIndex<2> &timeIdx) {
    calcMeridionalWindAdvection(timeIdx);
    calcMeridionalWindCoriolis(timeIdx);
    calcMeridionalWindPressureGradient(timeIdx);
} // calcMeridionalWindTendency

/**
 *  Input: u, v, ut
 *  Output, s1
 */
void BarotropicModel_A_ImplicitMidpoint::
calcZonalWindAdvection(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
            fu(i, j) = ut(timeIdx, i, j)*u(timeIdx, i, j);
            fv(i, j) = ut(timeIdx, i, j)*v(timeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            double dx1 = fu(i+1, j)-fu(i-1, j);
            double dy1 = fv(i, j+1)-fv(i, j-1);
            double dx2 = u(timeIdx, i, j)*(ut(timeIdx, i+1, j)-ut(timeIdx, i-1, j));
            double dy2 = v(timeIdx, i, j)*(ut(timeIdx, i, j+1)-ut(timeIdx, i, j-1))*cosLat[j];
            dut(i, j) = 0.5*((dx1+dx2)*factorLon[j]+(dy1+dy2)*factorLat[j]);
        }
    }
} // calcZonalWindAdvection

/**
 *  Input: u, v, vt
 *  Output, dvt
 */
void BarotropicModel_A_ImplicitMidpoint::
calcMeridionalWindAdvection(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL)-1; i <= mesh().ie(FULL)+1; ++i) {
            fu(i, j) = vt(timeIdx, i, j)*u(timeIdx, i, j);
            fv(i, j) = vt(timeIdx, i, j)*v(timeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            double dx1 = fu(i+1,j)-fu(i-1,j);
            double dy1 = fv(i,j+1)-fv(i,j-1);
            double dx2 = u(timeIdx, i, j)*(vt(timeIdx, i+1, j)-vt(timeIdx, i-1, j));
            double dy2 = v(timeIdx, i, j)*(vt(timeIdx, i, j+1)-vt(timeIdx, i, j-1))*cosLat[j];
            dvt(i, j) = 0.5*((dx1+dx2)*factorLon[j]+(dy1+dy2)*factorLat[j]);
        }
    }
} // calcMeridionalWindAdvection

/**
 *  Input: u, vt
 *  Output: dut
 */
void BarotropicModel_A_ImplicitMidpoint::
calcZonalWindCoriolis(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            double f = factorCor[j]+u(timeIdx, i, j)*factorCur[j];
            dut(i, j) -= f*vt(timeIdx, i, j);
        }
    }
} // calcZonalWindCoriolis

/**
 *  Input: u, ut
 *  Output: dvt
 */
void BarotropicModel_A_ImplicitMidpoint::
calcMeridionalWindCoriolis(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            double f = factorCor[j]+u(timeIdx, i, j)*factorCur[j];
            dvt(i, j) += f*ut(timeIdx, i, j);
        }
    }
} // calcMeridionalWindCoriolis

/*
 *  Input: gd, ghs, gdt
 *  Output: dut
 */
void BarotropicModel_A_ImplicitMidpoint::
calcZonalWindPressureGradient(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            dut(i, j) += (gd(timeIdx, i+1, j)-gd(timeIdx, i-1, j)+
                          ghs(i+1, j)-ghs(i-1, j))*
                         factorLon[j]*gdt(timeIdx, i, j);
        }
    }
} // calcZonalWindPressureGradient

/*
 *  Input: gd, ghs, gdt
 *  Output: dvt
 */
void BarotropicModel_A_ImplicitMidpoint::
calcMeridionalWindPressureGradient(const TimeLevelIndex<2> &timeIdx) {
    for (int j = mesh().js(FULL)+1; j <= mesh().je(FULL)-1; ++j) {
        for (int i = mesh().is(FULL); i <= mesh().ie(FULL); ++i) {
            dvt(i, j) += (gd(timeIdx, i, j+1)-gd(timeIdx, i, j-1)+
                          ghs(i, j+1)-ghs(i, j-1))*
                         factorLat[j]*cosLat[j]*gdt(timeIdx, i, j);
        }
    }
} // calcMeridionalWindPressureGradient

} // barotropic_model

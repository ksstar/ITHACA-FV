/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


/// \file
/// Source file of the unsteadyNS class.

#include "unsteadyNS_pimple.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
unsteadyNS_pimple::unsteadyNS_pimple() {}

unsteadyNS_pimple::unsteadyNS_pimple(int argc, char* argv[])
    :
    unsteadyNS(argc, argv)
{
    Info << offline << endl;
}

fvVectorMatrix unsteadyNS_pimple::get_Umatrix(volVectorField& U,
        volScalarField& p)
{
#include "initContinuityErrs.H"
    //Time& runTime = _runTime();
    IOMRFZoneList& MRF = _MRF();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix Ueqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
        + MRF.DDt(U)
        + turbulence->divDevReff(U)
        ==
        fvOptions(U)
    );
    Ueqn.relax();
    fvOptions.constrain(Ueqn);
    Ueqn_global = &Ueqn;
    return Ueqn;
}

fvScalarMatrix unsteadyNS_pimple::get_Pmatrix(volVectorField& U,
        volScalarField& p)
{
    IOMRFZoneList& MRF = _MRF();
    surfaceScalarField& phi = _phi();
    pimpleControl& pimple = _pimple();
    fv::options& fvOptions = _fvOptions();
    MRF.correctBoundaryVelocity(U);
    fvMesh& mesh = _mesh();
    volScalarField rAU(1.0 / Ueqn_global->A());
    volVectorField HbyA(constrainHbyA(rAU * Ueqn_global->H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA)+  		  fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);
    tmp<volScalarField> rAtU(rAU);

    if (pimple.consistent())
    {
        rAtU = 1.0 / max(1.0 / rAU - Ueqn_global->H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf();
        HbyA -= (rAU - rAtU()) * fvc::grad(p);
    }

    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);
	pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter()))); 
        
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    p.relax();
    U = HbyA - rAtU() * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );
    return pEqn;
}

void unsteadyNS_pimple::truthSolve2(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    IOMRFZoneList& MRF = _MRF();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // save initial condition in folder 0
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U);
    Pfield.append(p);
    counter++;
    nextWrite += writeEvery;


   // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
		fvVectorMatrix UEqn(get_Umatrix(U, p));

	    if (pimple.momentumPredictor())
	    {
    		solve(UEqn == -fvc::grad(p));
    		fvOptions.correct(U);
	    }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                get_Pmatrix(U, p);
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (int i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }

        runTime++;
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}

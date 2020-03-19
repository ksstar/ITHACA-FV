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
/// Source file of the unsteadyNSPISO class.


#include "unsteadyNSPISO.H"


unsteadyNSPISO::unsteadyNSPISO() {}

// Construct from zero
unsteadyNSPISO::unsteadyNSPISO(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
#include "fvOptions.H"
    _piso = autoPtr<pisoControl>
            (
                new pisoControl
                (
                    mesh
                )
            );
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
#include "createFields.H"
#include "createFvOptions.H"
    para = new ITHACAparameters;
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none", 
             "The BC method must be set to lift or penalty in ITHACAdict");
    timedepbcMethod = ITHACAdict->lookupOrDefault<word>("timedepbcMethod", "no");
    M_Assert(timedepbcMethod == "yes" || timedepbcMethod == "no",
             "The BC method can be set to yes or no");
    timeDerivativeSchemeOrder =
        ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");
    M_Assert(timeDerivativeSchemeOrder == "first"
             || timeDerivativeSchemeOrder == "second",
             "The time derivative approximation must be set to either first or second order scheme in ITHACAdict");
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
#include "fvCFD.H"

void unsteadyNSPISO::truthSolve(List<scalar> mu_now, fileName folder)

{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
    fv::options& fvOptions = _fvOptions();
#include "initContinuityErrs.H"
    pisoControl& piso = _piso();
    volScalarField& p = _p();
    volVectorField& U = _U();
    dimensionedScalar nu = _nu();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;

    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(phi, name(counter), folder);
    phi = linearInterpolate(U) & mesh.Sf();
    volScalarField Udiv = fvc::div(U);
    volScalarField Phidiv = fvc::div(phi);
    ITHACAstream::exportSolution(phi, name(counter), folder);
    ITHACAstream::exportSolution(Udiv, name(counter), folder);
    ITHACAstream::exportSolution(Phidiv, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U);
    Pfield.append(p);
    Phifield.append(phi);
    Udivfield.append(Udiv);
    Phidivfield.append(Phidiv);


    counter++;
    nextWrite += writeEvery;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- Pressure-velocity PISO corrector loop
        while (piso.correct())
        {
            {
#include "pEqn.H"
            }
        }


        //runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);

	    phi = linearInterpolate(U) & mesh.Sf();
	    volScalarField Udiv = fvc::div(U);
	    volScalarField Phidiv = fvc::div(phi);
	    ITHACAstream::exportSolution(phi, name(counter), folder);
	    ITHACAstream::exportSolution(Udiv, name(counter), folder);
	    ITHACAstream::exportSolution(Phidiv, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
	    Phifield.append(phi);
	    Udivfield.append(Udiv);
	    Phidivfield.append(Phidiv);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);

	}

    }

    fvVectorMatrix UEqn
    (
	fvm::ddt(U)
 	+ fvm::div(phi, U)
	- fvm::laplacian(nu, U)
     );
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);
    A = A;
    b = b;
   // Eigen::MatrixXd bw;
    bw.resize(b.size(),1);
    bw.col(0) = b;

     ITHACAstream::exportMatrix(bw, "bw", "eigen",
                               "./ITHACAoutput/BCvector/");

 ITHACAstream::SaveDenseMatrix(bw, "./ITHACAoutput/BCvector/",
                                  "bw");
 	    
    DivNorm.resize(Udivfield.size(),4);
    
    for (label j = 0; j < Udivfield.size(); j++)
    {
	DivNorm(j,0) = ITHACAutilities::L2norm(Udivfield[j]); //scalar div(U)
	DivNorm(j,1) = ITHACAutilities::L2norm(Phidivfield[j]); //scalar div(phi)
   	DivNorm(j,2) = runTime.deltaTValue()*
        	mag(Phidivfield[j])().weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,3) = runTime.deltaTValue()*
        	Phidivfield[j].weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar globalContErr
    }

    ITHACAstream::exportMatrix(DivNorm, "DivNorm", "eigen",
                                   folder);

}

bool unsteadyNSPISO::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}







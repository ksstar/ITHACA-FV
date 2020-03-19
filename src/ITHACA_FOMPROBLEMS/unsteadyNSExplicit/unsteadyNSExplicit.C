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
/// Source file of the unsteadyNSExplicit class.


#include "unsteadyNSExplicit.H"


unsteadyNSExplicit::unsteadyNSExplicit() {}

// Construct from zero
unsteadyNSExplicit::unsteadyNSExplicit(int argc, char* argv[])
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

void unsteadyNSExplicit::truthSolve(List<scalar> mu_now, fileName folder)

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
    dimensionedScalar dt = _dt();
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
    ITHACAstream::exportSolution(phi, name(counter), folder);

    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U);
    Pfield.append(p);
    Phifield.append(phi);

    counter++;
    nextWrite += writeEvery;

 /*   volVectorField  dU("dU", U * 0);
 volVectorField  U1("U1", U * 0);
 volVectorField  Un("Un", U * 0);
    volVectorField  U2("U2", U * 0);
    volVectorField  U3("U3", U * 0);
    volVectorField  U4("U4", U * 0);

    volScalarField  p2("p2", p * 0);
    volScalarField  p3("p3", p * 0);
    volScalarField  p4("p4", p * 0);


    surfaceScalarField phi1("phi1", phi*0);
    surfaceScalarField phi2("phi2", phi*0);
    surfaceScalarField phi3("phi3", phi*0);
    surfaceScalarField phi4("phi4", phi*0);*/

//

#include "CreatePoissonMatrix.H"
    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;
        
	phi = linearInterpolate(U) & mesh.Sf();
	volVectorField& Un = U; volVectorField& U1 = Un; surfaceScalarField& phi1 = phi; //Un = U; U1 = Un; phi1 = phi;
	
	// Stage 2
	volVectorField dU =      dt * a21 * (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	U = U1 + dU; 
	U.correctBoundaryConditions();
	phi = linearInterpolate(U) & mesh.Sf();
	solve(pEqn ==  (1/(c2*dt))*fvc::div(phi));

	U = U - (c2*dt)*fvc::grad(p);
	U.correctBoundaryConditions();
 	phi = linearInterpolate(U) & mesh.Sf();
	volVectorField& U2 = U; surfaceScalarField& phi2 = phi;

	//Stage 3
	dU =       dt * a31 *(-fvc::div(phi1,U1) + nu *fvc::laplacian(U1))
		 + dt * a32 *(-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	U = U1 +  dU;
	U.correctBoundaryConditions();
	phi = linearInterpolate(U) & mesh.Sf();
	solve(pEqn ==  (1/(c3*dt))*fvc::div(phi));

	U = Un - (c3*dt)*fvc::grad(p);
	U.correctBoundaryConditions();
 	phi = linearInterpolate(U) & mesh.Sf();
	volVectorField& U3 = U; surfaceScalarField& phi3 = phi;

	// Stage New Time Step
	dU =      dt * b1 *(-fvc::div(phi1,U1) + nu *fvc::laplacian(U1))
		+ dt * b2 *(-fvc::div(phi2,U2) + nu *fvc::laplacian(U2))
		+ dt * b3 *(-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	U = U1 +  dU; 
	U.correctBoundaryConditions();
	phi = linearInterpolate(U) & mesh.Sf();
	solve(pEqn ==  (1/dt)*fvc::div(phi));

	U = U - dt*fvc::grad(p);
	U.correctBoundaryConditions();
 	phi = linearInterpolate(U) & mesh.Sf();

	//U.storeOldTime();
	// U.store();

	// Stage pressure postProcessing

	//solve(pEqn ==  fvc::div(-fvc::div(phi,U) + nu *fvc::laplacian(U)));

	solve(pEqn ==  fvc::div(-fvc::div(phi,U)) );
	
	#include "continuityErrs.H"
        // Calculate the flow-direction filter tensor
    //    volScalarField magSqrU(magSqr(U));
    //    volSymmTensorField F(sqr(U)/(magSqrU + SMALL*average(magSqrU)));

        // Calculate the divergence of the flow-direction filtered div(U*U)
        // Filtering with the flow-direction generates a more reasonable
        // pressure distribution in regions of high velocity gradient in the
        // direction of the flow
      /*  volScalarField divDivUU
        (
            fvc::div
            (
                F & (-fvc::div(phi, U) + nu *fvc::laplacian(U)),
                "div(div(phi,U))"
            )
        );


	fvScalarMatrix pEqn
	(
	fvm::laplacian(p) - divDivUU
	);
	pEqn.setReference(pRefCell, pRefValue);

	pEqn.solve();*/

	


/*	volVectorField Uold = U; volVectorField Uc = U;
	phi = linearInterpolate(U) & mesh.Sf();

	// Stage 2
	dU =  dt* a21 *(-fvc::div(phi,U) + nu *fvc::laplacian(U));
	Uc = U + dU; U = Uold +0.5*dU;
    	#include "PressureCorrection.H"
 	phi = linearInterpolate(U) & mesh.Sf();

	//Stage 3
	dU =  dt * a31 *(-fvc::div(phi,U) + nu *fvc::laplacian(U));
		//dt  * a32 *(-fvc::div(phi,U2) + nu *fvc::laplacian(U2));
	Uc = U +  dU; U = Uold +0.5*dU;
    	#include "PressureCorrection.H"
 	phi = linearInterpolate(U) & mesh.Sf();

	//Stage 4
	dU =    dt * b1 *(-fvc::div(phi,U) + nu *fvc::laplacian(U));
		//dt * b2 *(-fvc::div(phi,U2) + nu *fvc::laplacian(U2))+
		//dt * b3 *(-fvc::div(phi,U3) + nu *fvc::laplacian(U3))	;

	Uc = U +  dU; U = Uold +0.5*dU;
    	#include "PressureCorrection.H"
 	phi = linearInterpolate(U) & mesh.Sf();/*

/*	 U2.correctBoundaryConditions();
	phi2 = linearInterpolate(U2) & mesh.Sf();
	solve(fvm::laplacian(p) == fvc::div(phi2));
	U2 = U2 - fvc::grad(p);
        U2.correctBoundaryConditions();


	//Stage 3

	dU =  dt * a31 *(-fvc::div(phi,U) + nu *fvc::laplacian(U))+
		dt  * a32 *(-fvc::div(phi,U2) + nu *fvc::laplacian(U2));
	U3 = U +  dU;
	 U3.correctBoundaryConditions();
	phi3 = linearInterpolate(U3) & mesh.Sf();
	solve(fvm::laplacian(p) == fvc::div(phi3));
	U3 = U3 - fvc::grad(p);
        U3.correctBoundaryConditions();


	//Stage 4

	dU =    dt * b1 *(-fvc::div(phi,U) + nu *fvc::laplacian(U))+
		dt  * b2 *(-fvc::div(phi,U2) + nu *fvc::laplacian(U2))+
		dt * b3 *(-fvc::div(phi,U3) + nu *fvc::laplacian(U3))	;

	U4 = U +  dU;
	 U4.correctBoundaryConditions();
	phi4 = linearInterpolate(U4) & mesh.Sf();
	//p = p;
	solve(fvm::laplacian(p) == fvc::div(phi4));
	U4 = U4 - fvc::grad(p);
        U4.correctBoundaryConditions();
	U = U4;*/



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

	volTensorField gradUt = fvc::grad(U);
//volScalarField divUt = gradUt.component(0) + gradUt.component(4) + gradUt.component(8); (Div is the trace of Grad, so...)

	   // volScalarField Phidiv = fvc::div(phi);

volScalarField Phidiv 
        (
		"Phidiv",
                (gradUt.component(0) + gradUt.component(4) + gradUt.component(8))
                

        );

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

DivNorm.resize(Udivfield.size(),2);
    
    for (label j = 0; j < Udivfield.size(); j++)
    {

   	DivNorm(j,0) = runTime.deltaTValue()*
        	mag(Phidivfield[j])().weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,1) = runTime.deltaTValue()*
        	Phidivfield[j].weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar globalContErr
    }

    ITHACAstream::exportMatrix(DivNorm, "DivNorm", "eigen",
                                   folder);

    

}

bool unsteadyNSExplicit::checkWrite(Time& timeObject)
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







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
    Method = ITHACAdict->lookupOrDefault<word>("Method", "RK3");
    M_Assert(Method == "FE" || Method == "RK3" , 
             "The BC method must be set to FE or RK3 in ITHACAdict");
    ExplicitMethod = ITHACAdict->lookupOrDefault<word>("ExplicitMethod", "A");
    M_Assert(ExplicitMethod == "Ales" || ExplicitMethod == "A" || ExplicitMethod == "B", 
             "The Explicit method must be set to Ales or A or B in ITHACAdict");
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
    fvMesh& mesh = _mesh();
    volScalarField p = _p();
    volVectorField U = _U();
    surfaceScalarField phi = _phi() ;
#include "initContinuityErrs.H"
    dimensionedScalar nu = _nu();
    dimensionedScalar dt = timeStep*_dt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Determine bc vector
    fvVectorMatrix UEqn
    (
	fvm::laplacian(nu, U)
    );     
    Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);
    bw.resize(b.size(),1);
    bw.col(0) = b;
    ITHACAstream::exportMatrix(bw, "bw0", "eigen",
                               "./ITHACAoutput/BCvector/");
    P0field.append(p);
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(phi, name(counter), folder);

    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U);
    Pfield.append(p);
    Phifield.append(phi);

    counter++;
    nextWrite += writeEvery;


if (ExplicitMethod == "Ales")
{
    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        
	// Stage 1
	volVectorField U1 = U; surfaceScalarField phi1 = 
	fvc::interpolate(U)& mesh.Sf(); volScalarField p1("p1", p*0);

	if (Method == "RK3")
	{
	    // Stage 2
	    p = p1;
	    volVectorField F1 = (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	    volVectorField U2 = U1 + dt * a21 *F1 ; 
	    U2.correctBoundaryConditions();

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c2*dt))*fvc::div(U2)); 
	    U2 = U2 - (c2*dt)*fvc::grad(p);
	    U2.correctBoundaryConditions();
	    surfaceScalarField phi2 = fvc::interpolate(U2)& mesh.Sf(); 
	    volScalarField p2("p2", p);
 	
	    //Stage 3
	    p = p1;
	    volVectorField F2  = (-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	    volVectorField U3 = U1 + dt * a31 * F1 + dt * a32 * F2  ;
	    U3.correctBoundaryConditions();

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c3*dt))*fvc::div(U3));
	    U3 = U3 - (c3*dt)*fvc::grad(p);
	    U3.correctBoundaryConditions();
	    surfaceScalarField phi3 = fvc::interpolate(U3)& mesh.Sf(); 
	    volScalarField p3("p3", p);
 	
	    // Stage New Time Step
	    p = p1;
	    volVectorField F3 =  (-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	    U = U1 +  dt * b1 * F1 +  dt * b2 * F2 + dt * b3 * F3 ; 
	    U.correctBoundaryConditions();

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U));
	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();
	}
	else if (Method == "FE")
	{
	    volVectorField F =  (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1) );
	    U = U1 +  dt *F  ; 
	    U.correctBoundaryConditions();

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U));
	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();

	}

	#include "continuityErrs.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
	    volVectorField gradP = fvc::grad(p);
            ITHACAstream::exportSolution(U, name(counter), folder);
	    ITHACAstream::exportSolution(gradP, name(counter), folder);
	    
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
	    Phifield.append(phi);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);

	}

    }
}
else if (ExplicitMethod == "A")
{

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        
	// Stage 1
	volVectorField U1 = U; surfaceScalarField phi1 = phi; 
	volScalarField p1("p1", p*0);

	if (Method == "RK3")
	{
	    // Stage 2
	    p = p1;
	    volVectorField F1 = (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	    volVectorField U2 = U1 + dt * a21 *F1 ; 
	    U2.correctBoundaryConditions();
	
	    surfaceScalarField phi2 = fvc::interpolate(U2)& mesh.Sf(); 

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c2*dt))*fvc::div(phi2)); 
	    U2 = U2 - (c2*dt)*fvc::grad(p);
	    U2.correctBoundaryConditions();

	    phi2 = phi2 - c2*dt*fvc::snGrad(p)*mag(mesh.magSf());
	    volScalarField p2("p2", p);
 	
	    //Stage 3
	    p = p1;
	    volVectorField F2  = (-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	    volVectorField U3 = U1 + dt * a31 * F1 + dt * a32 * F2  ;
	    U3.correctBoundaryConditions();

	    surfaceScalarField phi3 = fvc::interpolate(U3)& mesh.Sf(); 

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c3*dt))*fvc::div(phi3));
	    U3 = U3 - (c3*dt)*fvc::grad(p);
	    U3.correctBoundaryConditions();
	
	    phi3 = phi3 - c3*dt*fvc::snGrad(p)*mag(mesh.magSf());

	    volScalarField p3("p3", p);
 	
	    // Stage New Time Step
	    p = p1;
	    volVectorField F3 =  (-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	    U = U1 +  dt * b1 * F1 +  dt * b2 * F2 + dt * b3 * F3 ; 
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(phi));

	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();

	    phi = phi -  dt*fvc::snGrad(p)*mag(mesh.magSf());
	}
	else if (Method == "FE")
	{
	    p = p1;
	    volVectorField F =  (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1) );
	    U = U1 +  dt *F  ; 
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(phi));
	    p.correctBoundaryConditions();
	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();
	    phi = phi - dt*fvc::snGrad(p)*mag(mesh.magSf());

	}

	#include "continuityErrs.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
	    ITHACAstream::exportSolution(p, name(counter), folder);
	    ITHACAstream::exportSolution(phi, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
	    Phifield.append(phi);
	    //Phidiv = fvc::div(phi);
	    //Phidivfield.append(Phidiv);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);

	}

    }

    /*DivNorm.resize(Phidivfield.size(),2);
    for (label j = 0; j < Phidivfield.size(); j++)
    {
   	DivNorm(j,0) = runTime.deltaTValue()*
        	mag(Phidivfield[j])().weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,1) = runTime.deltaTValue()*
        	Phidivfield[j].weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar globalContErr
    }
    ITHACAstream::exportMatrix(DivNorm, "DivNorm", "eigen",folder);*/
}


else if (ExplicitMethod == "B")
{

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        
	// Stage 1
	volVectorField U1 = U; surfaceScalarField phi1 = phi; 
	volScalarField p1("p1", p*0);

	if (Method == "RK3")
	{
	    // Stage 2
	    p = p1;
	    volVectorField F1 = (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	    volVectorField U2 = U1 + dt * a21 *F1 ; 
	    U2.correctBoundaryConditions();
	
	    surfaceScalarField phi2 = phi1 + dt * a21 *fvc::flux(F1); 

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c2*dt))*fvc::div(phi2)); 
	    U2 = U2 - (c2*dt)*fvc::grad(p);
	    U2.correctBoundaryConditions();

	    //phi2 = phi2 - (c2)*( dt*fvc::snGrad(p)*mag(mesh.magSf()));
	    phi2 = phi2 - c2*dt*pEqn.flux();
	    volScalarField p2("p2", p);
 	
	    //Stage 3
	    p = p1;
	    volVectorField F2  = (-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	    volVectorField U3 = U1 + dt * a31 * F1 + dt * a32 * F2  ;
	    U3.correctBoundaryConditions();

	    surfaceScalarField phi3 = phi1 + dt * a31 * fvc::flux(F1) +
					dt * a32 * fvc::flux(F2); 

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c3*dt))*fvc::div(phi3));
	    U3 = U3 - (c3*dt)*fvc::grad(p);
	    U3.correctBoundaryConditions();
	
	    //phi3 = phi3 - (c3)*( dt*fvc::snGrad(p)*mag(mesh.magSf()));
	    phi3 = phi3 - c3*dt*pEqn.flux();
	    volScalarField p3("p3", p);
 	
	    // Stage New Time Step
	    p = p1;
	    volVectorField F3 =  (-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	    U = U1 +  dt * b1 * F1 +  dt * b2 * F2 + dt * b3 * F3 ; 
	    U.correctBoundaryConditions();

	    phi = phi1 + dt * b1 * fvc::flux(F1) +  dt * b2 * fvc::flux(F2)+ 
			dt * b3 * fvc::flux(F3); 

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(phi));

	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();

	    //phi = phi - ( dt*fvc::snGrad(p)*mag(mesh.magSf()));
	    phi = phi - dt*pEqn.flux();
	}
	else if (Method == "FE")
	{
	    p = p1;
	    volVectorField F =  (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1) );
	    U = U1 +  dt *F  ; 
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(phi));
	    p.correctBoundaryConditions();
	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();
	    //phi = phi - ( dt*fvc::snGrad(p)*mag(mesh.magSf()));
 	    phi = phi - dt*pEqn.flux();
	}

	#include "continuityErrs.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
	    ITHACAstream::exportSolution(p, name(counter), folder);
	    ITHACAstream::exportSolution(phi, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
	    Phifield.append(phi);
	    //Phidiv = fvc::div(phi);
	    //Phidivfield.append(Phidiv);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);

	}

    }

/*    DivNorm.resize(Phidivfield.size(),2);
    for (label j = 0; j < Phidivfield.size(); j++)
    {
   	DivNorm(j,0) = runTime.deltaTValue()*
        	mag(Phidivfield[j])().weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,1) = runTime.deltaTValue()*
        	Phidivfield[j].weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar globalContErr
    }
    ITHACAstream::exportMatrix(DivNorm, "DivNorm", "eigen",
                                   folder);*/
}


   
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

   





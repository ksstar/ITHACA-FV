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
    M_Assert(Method == "FE" || Method == "RK3" || Method == "midPoint", 
             "The Method must be set to FE or RK3 or midPoint in ITHACAdict");
    ExplicitMethod = ITHACAdict->lookupOrDefault<word>("ExplicitMethod", "A");
    M_Assert(ExplicitMethod == "Ales" || ExplicitMethod == "A" || ExplicitMethod == "B", 
             "The Explicit method must be set to Ales or A or B in ITHACAdict");
    PoissonMethod = ITHACAdict->lookupOrDefault<word>("PoissonMethod", "A");
    M_Assert(PoissonMethod == "FOM" || PoissonMethod == "ROM", 
             "The Poisson method must be set to FOM or ROM in ITHACAdict");
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
    Vector<double> inl(0, 0, 0);
    volVectorField U0 = _U();
    assignIF(U0, inl);
   ITHACAutilities::changeBCtype( U0,"fixedValue",1);
    	    assignBC( U0,1,inl);

dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(0.005)
        );

dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(0.01)
        );

   /* volScalarField divU = (1/dt_fake)*fvc::div(U0);
    ITHACAstream::exportSolution(divU, name(counter), folder);
    

    volScalarField divlap = fvc::div(fvc::laplacian(nu_fake, U0));
    ITHACAstream::exportSolution(divlap, name(counter), folder);*/

    surfaceScalarField phi("phi",_phi()) ;
    //surfaceScalarField phi0("phi0",_phi()) ;
    if (ExplicitMethod == "Ales")
    {
    //    phi0 = fvc::flux(U0);//_phi();
	phi = fvc::flux(U);
    }
    else if (ExplicitMethod == "A")
    {
       // phi0 = _phi();
	phi = _phi();
    }
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
	fvm::laplacian(nu, U0)
	-fvm::div(fvc::flux(U0),U0)
    );     

    Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);
    bw.resize(b.size(),1);
    bw.col(0) = b;
    ITHACAstream::exportMatrix(bw, "bw0", "eigen",
                               "./ITHACAoutput/BCvector/");
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
	volVectorField U1(U.oldTime()); surfaceScalarField phi1 = 
	fvc::interpolate(U1)& mesh.Sf(); volScalarField p1(p);

	if (Method == "RK3")
	{
	    // Stage 2
	    p = p*0;
	    p.correctBoundaryConditions();
	    volVectorField F1aux(U);
	    F1aux =  dt * (a21 *(fvc::laplacian(nu,U1)-fvc::div(phi1,U1)));

	    volVectorField U2(U);
	    U2 = U.oldTime() + F1aux;
	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/(c2*dt))*fvc::div(U2) ); 

	    U2 = U2 - c2*dt*fvc::grad(p);

	    surfaceScalarField phi2 = fvc::interpolate(U2)& mesh.Sf(); 
 	
	    //Stage 3
	    p = p*0;
	    p.correctBoundaryConditions();	
	    volVectorField F2aux(U);
	    F2aux =  dt * (a31 * (fvc::laplacian(nu,U1)-fvc::div(phi1,U1)) +
			   a32 * (fvc::laplacian(nu,U2)-fvc::div(phi2,U2)));

	    volVectorField U3(U);
            U3 = U.oldTime() + F2aux ;
	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/(c3*dt))*fvc::div(U3));

	    U3 = U3 - c3*dt*fvc::grad(p);
	    
	    surfaceScalarField phi3 = fvc::interpolate(U3)& mesh.Sf(); 
 	
	  /*  // Stage 4
	    p = p*0;
	    p.correctBoundaryConditions();	
	    volVectorField F4aux(U);
	    F4aux =  dt * (a41 * (fvc::laplacian(nu,U1)-fvc::div(phi1,U1)) +
			   a42 * (fvc::laplacian(nu,U2)-fvc::div(phi2,U2)) +
			   a43 * (fvc::laplacian(nu,U3)-fvc::div(phi3,U3))
			  );
 	
	    volVectorField U4(U);
	    U4 = U.oldTime() + F4aux;
	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/(c4*dt))*fvc::div(U4) );
	    U4 = U4 - c4*dt*fvc::grad(p);

	    surfaceScalarField phi4 = fvc::interpolate(U4)& mesh.Sf();*/

	    // Stage New Time Step
	    p = p*0;
	    p.correctBoundaryConditions();	
	    volVectorField F3aux(U);
	    F3aux =  dt * (b1 * (fvc::laplacian(nu,U1)-fvc::div(phi1,U1)) +
			   b2 * (fvc::laplacian(nu,U2)-fvc::div(phi2,U2)) +
			   b3 * (fvc::laplacian(nu,U3)-fvc::div(phi3,U3)) //+
			   //b4 * (fvc::laplacian(nu,U4)-fvc::div(phi4,U4))
			  );
 	
	    volVectorField Uaux(U);
	    Uaux = U.oldTime() + F3aux;
	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/dt)*fvc::div(Uaux) );
	    U = Uaux - dt*fvc::grad(p);
	    Ufield.append(U);
	    phi = fvc::interpolate(U)& mesh.Sf();

	}
	else if (Method == "midPoint")
	{
	    // Stage 2
	    volVectorField F1aux(U);
	    F1aux =  dt * (a21 *(fvc::laplacian(nu,U1)));

	    volVectorField U2(U);
	    U2 = U.oldTime() + F1aux;
	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/(c2*dt))*fvc::div(U2) ); 

	    U2 = U2 - c2*dt*fvc::grad(p);

	    surfaceScalarField phi2 = fvc::interpolate(U2)& mesh.Sf(); 
 

	    // Stage New Time Step
	    p = p*0;
	    p.correctBoundaryConditions();	
	    volVectorField F3aux(U);
	    F3aux =  dt * (b1 * (fvc::laplacian(nu,U1)) +
			   b2 * (fvc::laplacian(nu,U2))
			  );
 	
	    volVectorField Uaux(U);
	    Uaux = U.oldTime() + F3aux;
	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1.0/dt)*fvc::div(Uaux) );
	    U = Uaux - dt*fvc::grad(p);

	    phi = fvc::interpolate(U)& mesh.Sf();

	    Ufield.append(U);

	}
	else if (Method == "FE")
	{
	    volVectorField U1aux(U);
	    volVectorField Faux(U);

	    Faux = dt * (fvc::laplacian(nu,U1)-fvc::div(phi,U1));
	    U1aux = U.oldTime() + Faux;

	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U1aux)); 

	    U = U1aux - dt*fvc::grad(p);

	    phi = fvc::interpolate(U)& mesh.Sf();
	    

	}

	#include "continuityErrs.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
	    ITHACAstream::exportSolution(p, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            
	    Ufield.append(U);
            Pfield.append(p);
	    Phifield.append(phi);
	
	    ITHACAstream::exportSolution(phi, name(counter), folder);

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
	volVectorField U1 = U.oldTime(); surfaceScalarField phi1=phi.oldTime(); 
	volScalarField p1(p);

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
	    //U2 = U2 - (c2*dt)*fvc::grad(p);
	    //U2.correctBoundaryConditions();

	    //phi2 = phi2 - c2*dt*fvc::snGrad(p)*mag(mesh.magSf());
	    //volScalarField p2("p2", p);
 	
	    //Stage 3
	    p = p1;
	    volVectorField F2  = (-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	    volVectorField U3 = U1 + dt * a31 * F1 + dt * a32 * F2  ;
	    U3.correctBoundaryConditions();

	    surfaceScalarField phi3 = fvc::interpolate(U3)& mesh.Sf(); 

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/(c3*dt))*fvc::div(phi3));
	   // U3 = U3 - (c3*dt)*fvc::grad(p);
	   // U3.correctBoundaryConditions();
	
	    //phi3 = phi3 - c3*dt*fvc::snGrad(p)*mag(mesh.magSf());

	    //volScalarField p3("p3", p);
 	
	    // Stage New Time Step
	    p = p1;
	    volVectorField F3 =  (-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	    U = U1 +  dt * b1 * F1 +  dt * b2 * F2 + dt * b3 * F3 ; 
	    U.correctBoundaryConditions();
	    UStarfield.append(U);
	    phi = fvc::interpolate(U)& mesh.Sf();

	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(phi));

	    //U = U - dt*fvc::grad(p);
	    //U.correctBoundaryConditions();

	   // phi = phi -  dt*fvc::snGrad(p)*mag(mesh.magSf());



	}
	else if (Method == "FE")
	{

	
	    volVectorField U1aux(U);
	    volVectorField Faux(U);

	    Faux = dt * (-fvc::div(phi1,U1)+ fvc::laplacian(nu,U1));
	    U1aux = U.oldTime() + Faux;
	    p = p*0;
	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U1aux)); 

	    U = U1aux - dt*fvc::grad(p);
	    U.correctBoundaryConditions();
	    phi = fvc::flux(U1aux) - dt*fvc::snGrad(p)*mag(mesh.magSf());
	    


	    /*volVectorField U1aux(U);
	    volVectorField Caux(U);
	    volVectorField Daux(U);
	    
	    U1aux = U1;
	    Daux = dt * nu *fvc::laplacian(U1);
	    Caux = dt * (-fvc::div(phi1,U1));

	    //phi = fvc::flux(U1aux+Daux+Caux);
	    fvScalarMatrix pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U1aux) +(1/dt)*fvc::div(Daux) +(1/dt)*fvc::div(Caux));
	    U = U1aux+Daux+Caux - dt*fvc::grad(p);
	    //
	    U.correctBoundaryConditions();
	    phi = fvc::flux(U1aux)+fvc::flux(Daux)+fvc::flux(Caux) - dt*fvc::snGrad(p)*mag(mesh.magSf());*/

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
	
	    surfaceScalarField phi2 = (phi1 + dt * a21 *fvc::flux(F1)); 

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
	    UStarfield.append(U);
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
	    UStarfield.append(U);
	    phi = phi1 + dt * fvc::flux(F);

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

   

/* Stage 2
	    p = p1;
	    volVectorField F1 = (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	    volVectorField U2 = U1 + dt * a21 *F1 ; volScalarField p1("p1", p*0);
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
	    UStarfield.append(U);
	    pEqn = fvm::laplacian(p);
	    pEqn.setReference(pRefCell, pRefValue);
	    solve(pEqn ==  (1/dt)*fvc::div(U));
	    U = U - dt*fvc::grad(p);
	    U.correctBoundaryConditions();

	    phi = fvc::interpolate(U)& mesh.Sf();*/



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
volScalarField pzero = p;
    volVectorField U = _U();
    surfaceScalarField phi = _phi();
#include "initContinuityErrs.H"
    dimensionedScalar nu = _nu();
    dimensionedScalar dt = _dt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Determine bc vector
  //  fvVectorMatrix UEqn
  //  (
//	fvm::laplacian(dimensionedScalar("1", dimless, 1), U)
  //  );     
  //  Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);
  //  bw.resize(b.size(),1);
  //  bw.col(0) = b;
 //   ITHACAstream::exportMatrix(bw, "bw0", "eigen",
 //                              "./ITHACAoutput/BCvector/");

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

	surfaceScalarField phiStarF1("phiStarF1", phi*0);
	surfaceScalarField phiStar("phiStar", phi*0);
	volScalarField phiF1 ("phiF1",fvc::div(phi)*0);
	volScalarField phiUminusU1  ("phiUminusU1",fvc::div(phi)*0);
	volScalarField DivphiStar("DivphiStar",fvc::div(phi)*0);
	volScalarField DivphiStarF1("DivphiStarF1",fvc::div(phi)*0);

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
	volVectorField U1 = U; surfaceScalarField phi1 = phi; volScalarField p1("p1", p);

	if (Method == "RK3")
	{
	// Stage 2
	p = p1;
	volVectorField F1 = (-fvc::div(phi1,U1) + nu *fvc::laplacian(U1));
	volVectorField U2 = U1 + dt * a21 *F1; 

	U2.correctBoundaryConditions();

	// Method A
	//surfaceScalarField phi2 = fvc::interpolate(U2)& mesh.Sf(); 

	// Method B1
	surfaceScalarField phiStar2 =   fvc::interpolate(U2-U1) & mesh.Sf();
	surfaceScalarField phi2 = phi1 +  phiStar2;

	/*// Method B2
	surfaceScalarField phiStar2 =   fvc::interpolate(dt * a21 *F1) & mesh.Sf() ;
	surfaceScalarField phi2 = phi1 + phiStar2;

	fvScalarMatrix pEqn = fvm::laplacian(phi2);
	pEqn.setReference(pRefCell, pRefValue);
	solve(pEqn ==  (1/(c2*dt))*fvc::div(phi2)); 
	p.correctBoundaryConditions();
	U2 = U2 - (c2*dt)*fvc::grad(p);
	U2.correctBoundaryConditions();

	phi2 = phi2 - (c2*dt)*pEqn.flux();
	volScalarField p2("p2", p);*/
 	
	//Stage 3
	p = p1;
	volVectorField F2  = (-fvc::div(phi2,U2) + nu *fvc::laplacian(U2));
	volVectorField U3 = U1 + dt * a31 * F1 + dt * a32 * F2  ;
	U3.correctBoundaryConditions();

	// Method A
	// surfaceScalarField phi3 = fvc::interpolate(U3)& mesh.Sf();

	// Method B1
	surfaceScalarField phiStar3 =   fvc::interpolate(U3-U1) & mesh.Sf() ;
	surfaceScalarField phi3 = phi1 + phiStar3;

/*	// Method B2
	surfaceScalarField phiStar3 =   fvc::interpolate(dt * a31 * F1 
				+ dt * a32 * F2) & mesh.Sf() ;
	surfaceScalarField phi3 = phi1 + phiStar3;

	pEqn = fvm::laplacian(p);
	pEqn.setReference(pRefCell, pRefValue);
	solve(pEqn ==  (1/(c3*dt))*fvc::div(phi3));
	p.correctBoundaryConditions();
	U3 = U3 - (c3*dt)*fvc::grad(p);
	U3.correctBoundaryConditions();

	phi3 = phi3 - (c3*dt)*pEqn.flux();
	volScalarField p3("p3", p);*/
 	
	// Stage New Time Step
	p = p1;
	volVectorField F3 =  (-fvc::div(phi3,U3) + nu *fvc::laplacian(U3));
	U = U1 +  dt * b1 * F1 +  dt * b2 * F2 + dt * b3 * F3 ; 
	U.correctBoundaryConditions();

	// Method A
	//phi = fvc::interpolate(U)& mesh.Sf();

	// Method B1
	surfaceScalarField phiStar =   fvc::interpolate(U-U1) & mesh.Sf() ;
	phi = phi1 + phiStar;	

/*	// Method B2
	//surfaceScalarField phiStar =   fvc::interpolate(dt * b1 * F2  
	//			+ dt * b2 * F3 + dt * b3 * F) & mesh.Sf() ;
	//phi = phi1 + phiStar;	

	pEqn = fvm::laplacian(p);
	pEqn.setReference(pRefCell, pRefValue);
	solve(pEqn ==  (1/dt)*fvc::div(phi));
	p.correctBoundaryConditions();
	U = U - dt*fvc::grad(p);
	U.correctBoundaryConditions();

        phi = phi - dt*pEqn.flux();*/
	}
	else if (Method == "FE")
	{

	p = p1;
	volVectorField Un("Un", U*0);
	volVectorField F1("F1", U*0);

	F1 =  dt * b1 *(-fvc::div(phi1,U1,"Gauss midPoint") + nu *fvc::laplacian(U1,"Gauss midPoint orthogonal"));
	Un = U1 + F1; 

	phiF1 = fvc::div(F1, "Gauss midPoint");	
	phiStarF1 =   fvc::interpolate(F1,"midPoint") & mesh.Sf() ;
	DivphiStarF1 =  fvc::div(phiStarF1) ;

	phiUminusU1 = fvc::div(Un, "Gauss midPoint")-fvc::div(U1, "Gauss midPoint");
	phiStar =   fvc::interpolate((Un-U1),"midPoint") & mesh.Sf() ;
	DivphiStar =  fvc::div(phiStar) ;

	// Method A
	//phi = fvc::interpolate(U)& mesh.Sf();

	// Method B1
	surfaceScalarField phin = phi1 + phiStar;	

	// Method B2
	//surfaceScalarField phiStar =   fvc::interpolate(F1) & mesh.Sf() ;
	
	

	fvScalarMatrix pEqn = fvm::laplacian(p);
	pEqn.setReference(pRefCell, pRefValue);
	solve(pEqn ==  ((1/dt)*(fvc::div(phin))));
	p.correctBoundaryConditions();
	U = Un - dt*fvc::grad(p);
	U.correctBoundaryConditions();

	phi = phin - dt * pEqn.flux();

	}

	#include "continuityErrs.H"


	

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        
 	if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
	    ITHACAstream::exportSolution(p, name(counter), folder);

	    //volScalarField Udiv = fvc::div(U);
	    //volScalarField Phidiv = fvc::div(phi);
	    ITHACAstream::exportSolution(phi, name(counter), folder);
	    ITHACAstream::exportSolution(phiStar, name(counter), folder);
	    ITHACAstream::exportSolution(phiStarF1, name(counter), folder);
	    ITHACAstream::exportSolution(phiF1, name(counter), folder);
	    ITHACAstream::exportSolution(phiUminusU1, name(counter), folder);
	    ITHACAstream::exportSolution(DivphiStar, name(counter), folder);
	    ITHACAstream::exportSolution(DivphiStarF1, name(counter), folder);
	    //ITHACAstream::exportSolution(Udiv, name(counter), folder);
	    //ITHACAstream::exportSolution(Phidiv, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
	    //Phifield.append(phi);
	   // Udivfield.append(Udiv);
	    //Phidivfield.append(Phidiv);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);

	}

    }

   /* DivNorm.resize(Udivfield.size(),4);
    for (label j = 0; j < Udivfield.size(); j++)
    {

   	DivNorm(j,0) = runTime.deltaTValue()*
        	mag(Phidivfield[j])().weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,1) = runTime.deltaTValue()*
        	Phidivfield[j].weightedAverage(Phidivfield[j].mesh().V()).value(); //scalar globalContErr

	DivNorm(j,2) = runTime.deltaTValue()*
        	mag(Udivfield[j])().weightedAverage(Udivfield[j].mesh().V()).value(); //scalar sumLocalContErr
	DivNorm(j,3) = runTime.deltaTValue()*
        	Udivfield[j].weightedAverage(Udivfield[j].mesh().V()).value(); //scalar globalContErr
    }
    ITHACAstream::exportMatrix(DivNorm, "DivNorm", "eigen",
                                   folder);*/

   
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







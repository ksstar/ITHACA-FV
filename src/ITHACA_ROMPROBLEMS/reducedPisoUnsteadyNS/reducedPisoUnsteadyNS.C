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
/// Source file of the reducedSteadyNS class

#include "reducedPisoUnsteadyNS.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedPisoUnsteadyNS::reducedPisoUnsteadyNS() {}

reducedPisoUnsteadyNS::reducedPisoUnsteadyNS(unsteadyNS_piso& FOMproblem)
    :
    problem(&FOMproblem)
{
    // Create a new Umodes set where the first one is the lift function
    ULmodes.append(problem->liftfield[0]);

    for (int i = 0; i < problem->Umodes.size(); i++)
    {
        ULmodes.append(problem->Umodes.toPtrList()[i]);
    }
}



// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedPisoUnsteadyNS::solveOnline_Piso()
{
    counter++;
    Eigen::VectorXd uresidualOld;
    Eigen::VectorXd presidualOld;
    uresidualOld.resize(ULmodes.size());
    presidualOld.resize(problem->Pmodes.size());
    Eigen::VectorXd uresidual;
    Eigen::VectorXd presidual;
    scalar residual_jump(1);
    scalar U_norm_res(1);
    scalar P_norm_res(1);
  
    Time& runTime = problem->_runTime();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime+dt);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(dt);

    label i = 1;

    Eigen::MatrixXd a = Eigen::VectorXd::Zero(ULmodes.size());
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(problem->Pmodes.size());

    volVectorField Uaux("Uaux", problem->Ufield[0]);
    Uaux.storeOldTime(); 
    volScalarField Paux("Paux", problem->Pfield[0]);
    Paux.storeOldTime(); 

    ITHACAstream::exportSolution(Uaux, name(i),
                                 "./ITHACAoutput/Reconstruct/");
    ITHACAstream::exportSolution(Paux, name(i),
                                 "./ITHACAoutput/Reconstruct/");
    i++;

    ITHACAparameters para;
    float residualJumpLim =
        para.ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
    float normalizedResidualLim =
        para.ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);

    pisoControl& piso = problem->_piso();
    setRefCell(Paux, piso.dict(), problem->pRefCell, problem->pRefValue);

    while (runTime.loop())
    {
	Info << "Time = " << runTime.timeName() << nl << endl;	

        residual_jump = 1;
        U_norm_res = 1;
        P_norm_res = 1;

 	 for (int i = 0; i < 5; i++)
	 {

		Uaux = ULmodes.reconstruct(a, "Uaux");
            	Paux = problem->Pmodes.reconstruct(b, "Paux");
            	fvVectorMatrix Au(get_Umatrix_Online(Uaux, Paux));
		
		List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(Au);
       		a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
        	Info << "res for a: " << uresidual.norm() << endl;
		
		uresidual = uresidual.cwiseAbs();
		U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
		//}

 		Uaux = ULmodes.reconstruct(a, "Uaux");
        	problem->_phi() = linearInterpolate(Uaux) & problem->_U().mesh().Sf();

		/// Construct pressure matrix using the momentum matrix
		while (piso.correct())
        	{
        		fvScalarMatrix Ap(get_Pmatrix_Online(Uaux, Paux));
			List<Eigen::MatrixXd> RedLinSysP = problem->Pmodes.project(Ap);
		b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
        	Info << "res for b: " << presidual.norm() << endl;
		}

 		Uaux = ULmodes.reconstruct(a, "Uaux");
        	Paux = problem->Pmodes.reconstruct(b, "Paux");
    	}
  
	Info << "Time = " << i*dt << nl << endl;
        ITHACAstream::exportSolution(Uaux, name(i),
                                 "./ITHACAoutput/Reconstruct/");
        ITHACAstream::exportSolution(Paux, name(i),
                                 "./ITHACAoutput/Reconstruct/");
	Uaux.storeOldTime(); 
        Paux.storeOldTime(); 
	
   	i++ ;  
    }
}

fvVectorMatrix reducedPisoUnsteadyNS::get_Umatrix_Online(volVectorField& U,
        volScalarField& p)
{
    Time& runTime = problem->_runTime();
    fvMesh& mesh = problem->_mesh();
    IOMRFZoneList& MRF = problem->_MRF();
    surfaceScalarField& phi = problem->_phi();
    MRF.correctBoundaryVelocity(U);

    IOdictionary transportProperties
    (
    	IOobject
    	(
        	"transportProperties",
        	runTime.constant(),
        	mesh,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
    	)
   );

dimensionedScalar nu
(
     transportProperties.lookup("nu")
);

    fvVectorMatrix Ueqn
    (
        fvm::ddt(U) 
	+ fvm::div(phi, U)
        -fvm::laplacian(nu,U) == -fvc::grad(p)
    );
    problem->Ueqn_global = &Ueqn;
    return Ueqn;
}

fvScalarMatrix reducedPisoUnsteadyNS::get_Pmatrix_Online(volVectorField& U,
        volScalarField& p)
{
    surfaceScalarField& phi = problem->_phi();
    fvMesh& mesh = problem->_mesh();
    
   pisoControl piso(mesh);
    
    volScalarField rAU(1.0 / problem->Ueqn_global->A());
    volVectorField HbyA(constrainHbyA(rAU * problem->Ueqn_global->H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA)+  		  fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));
    adjustPhi(phiHbyA, U, p);
// Update the pressure BCs to ensure flux consistency    
    constrainPressure(p, U, phiHbyA, rAU);
//while (piso.correctNonOrthogonal()) 
//    {  
//    fvScalarMatrix pEqn
//   (
//        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
//    );
//pEqn.setReference(problem->pRefCell, problem->pRefValue);
//pEqn.solve(mesh.solver(p.select(piso.finalInnerIter()))); 
//
//        if (piso.finalNonOrthogonalIter())
//        {
//            phi = phiHbyA - pEqn.flux();
//        }
//}

U = HbyA - rAU * fvc::grad(p);

    U.correctBoundaryConditions();

fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );

pEqn.setReference(problem->pRefCell, problem->pRefValue);

    return pEqn;
}


void reducedPisoUnsteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.size()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    vel_now = vel;
}



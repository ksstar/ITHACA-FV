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
Description
    Example of an unsteady NS Reduction Problem
SourceFiles
    19ExplicitOpen.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNSExplicit.H"
#include "ReducedUnsteadyNSExplicit.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial19ExplicitOpen: public unsteadyNSExplicit
{
    public:
        explicit tutorial19ExplicitOpen(int argc, char* argv[])
            :
            unsteadyNSExplicit(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");

		//ITHACAstream::read_fields(Ufield_on, U, "./ITHACAoutputError/OfflineAlesDiffRK3/");
		//ITHACAstream::read_fields(Pfield_on, p, "./ITHACAoutputError/OfflineAlesDiffRK3/");

 
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    //assignIF(U, inl);
                    change_viscosity( mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }


	void liftSolve()
{

    for (label k = 0; k < inletIndex.rows(); k++)
    {

        Time& runTime = _runTime();

	volScalarField p = _p();
        volVectorField U = _U();
	fvMesh& mesh = _mesh();
        surfaceScalarField phi= fvc::interpolate(U) & mesh.Sf();
        
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
	pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            
        }
	phi = fvc::interpolate(Ulift) & mesh.Sf();
        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);
        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();
            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(Ulift);
    }
}

};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
#include "fvCFD.H"
    // Construct the tutorial19ExplicitOpen object
    tutorial19ExplicitOpen example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 2.0;
    example.timeStep = 0.005;
    example.writeEvery = 0.005;

    if (example.Method == "RK3")
    {
    // RK coefficients (RK3)
    /*example.a21 = 0.33333333333333333333;
    example.a31 = -1.0;
    example.a32 = 2.0;
    example.b1 = 0;
    example.b2 = 0.75;
    example.b3 = 0.25;
    example.c1 = 0;
    example.c2 = 0.3333333333333333333333; 
    example.c3 = 1.0;*/

    example.a21 = 0.5;
    example.a31 = -1.0;
    example.a32 = 2.0;
    example.b1 = 0.1666666666666666666;
    example.b2 = 0.6666666666666666666;
    example.b3 = 0.1666666666666666666;
    example.c1 = 0;
    example.c2 = 0.5; 
    example.c3 = 1.0;

    /*example.a21 = 0.5;
    example.a31 = -1.0;
    example.a32 = 2.0;
    example.b1 = 1.0;
    example.b2 = 0.0;
    example.b3 = 0.0;
    example.c1 = 0;
    example.c2 = 0.5; 
    example.c3 = 1.0;*/

    /*example.a21 = 0.5;
    example.a31 = 0.5;
    example.a32 = 0.5;
    example.b1 = 0.3333333333333333333;
    example.b2 = 0.3333333333333333333;
    example.b3 = 0.3333333333333333333;
    example.c1 = 0;
    example.c2 = 0.5; 
    example.c3 = 1.0;*/

    /*example.a21 = 1.0;
    example.a31 = 0.375;
    example.a32 = 0.125;
    example.a41 = -0.125;
    example.a42 = -0.375;
    example.a43 = 1.5;
    example.b1 = 0.166666666666666666666666;
    example.b2 = -0.05555555555555555555555;
    example.b3 = 0.666666666666666666666666;
    example.b4 = 0.222222222222222222222222;
    example.c1 = 0;
    example.c2 = 1.0; 
    example.c3 = 0.5;
    example.c4 = 1.0;*/
    }
    else if (example.Method == "FE")
    {
    example.a21 = 0;
    example.a31 = 0;
    example.a32 = 0;
    example.b1 = 1.0;
    example.b2 = 0;
    example.b3 = 0;
    example.c1 = 0;
    example.c2 = 0;
    example.c3 = 0;
    }
    else if (example.Method == "midPoint")
    {
    example.a21 = 0.5;
    example.b1 = 0.0;
    example.b2 = 1.0;
    example.c1 = 0;
    example.c2 = 0.5;
    }

    // Perform The Offline Solve;
    auto start_FOM = std::chrono::high_resolution_clock::now();
    example.offlineSolve();
    auto finish_FOM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_FOM = finish_FOM - start_FOM;
//exit(0);

if (example.bcMethod == "lift")
    {
    // Search the lift function
     example.liftSolve();
    // Normalize the lifting function
  //  ITHACAutilities::normalizeFields(example.liftfield);
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
  ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    }
else
{    
	ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);

}



 /*  // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields_abs(
                                         example.Ufield_on, example.Ufield);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields_abs(
                                         example.Pfield_on, example.Pfield);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");

 // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixInfU = ITHACAutilities::error_listfields_inf(
                                         example.Ufield_on, example.Ufield);
    Eigen::MatrixXd L2errorMatrixInfP = ITHACAutilities::error_listfields_inf(
                                         example.Pfield_on, example.Pfield);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixInfU, "L2errorMatrixInfU", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");
    ITHACAstream::exportMatrix(L2errorMatrixInfP, "L2errorMatrixInfP", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");

 // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L1errorMatrixU = ITHACAutilities::error_listfields_Kazemi(
                                         example.Ufield_on, example.Ufield);
    Eigen::MatrixXd L1errorMatrixP = ITHACAutilities::error_listfields_Kazemi(
                                         example.Pfield_on, example.Pfield);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L1errorMatrixU, "L1errorMatrixU", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");
    ITHACAstream::exportMatrix(L1errorMatrixP, "L1errorMatrixP", "eigen",
                               "./ITHACAoutput/l2errorBdt0p1");
exit(0);*/


    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                       NmodesPout);

if (example.ExplicitMethod == "A"|| example.ExplicitMethod == "B")
{
    ITHACAPOD::getModesPhi(example.Phifield, example.Phimodes, example.podex, 0, 0,
                       NmodesUout);

    Eigen::MatrixXd  DivModes;
    DivModes.resize(example.Phimodes.size(),2);


    for (int j = 0; j < example.Phimodes.size(); j++)
    {

	adjustPhi(example.Phimodes[j], example.Ufield[0], example.Pfield[0]);
	volScalarField PhidivDivModes = fvc::div(example.Phimodes[j]);

	 

   	DivModes(j,0) = example.timeStep*
        	mag(PhidivDivModes)().weightedAverage(PhidivDivModes.mesh().V()).value(); //scalar sumLocalContErr
	DivModes(j,1) = example.timeStep*
        	PhidivDivModes.weightedAverage(PhidivDivModes.mesh().V()).value(); //scalar globalContErr


    }
    ITHACAstream::exportMatrix(DivModes, "DivModes", "eigen",
                                   "./ITHACAoutput/POD/");
}
//exit(0);
/*
    Eigen::MatrixXd  DivRec;
    DivRec.resize(example.Ufield.size(),2);

	Eigen::MatrixXd coeffL2U = ITHACAutilities::get_coeffs(example.Ufield,
              example.Umodes, 2);
	PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                example.Umodes, coeffL2U, 2);

    for (int j = 0; j < example.Ufield.size(); j++)
    {

	volScalarField UdivDivRec = fvc::div(rec_fieldU[j]);
	
   	DivRec(j,0) = example.timeStep*
        	mag(UdivDivRec)().weightedAverage(UdivDivRec.mesh().V()).value(); //scalar sumLocalContErr
	DivRec(j,1) = example.timeStep*
        	UdivDivRec.weightedAverage(UdivDivRec.mesh().V()).value(); //scalar globalContErr

    }
    ITHACAstream::exportMatrix(DivRec, "DivRec", "eigen",
                                   "./ITHACAoutput/POD/");*/


 // Perform the projection for all number of modes in list List_of_modes
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size(),1);
    Eigen::MatrixXd L2errorProjMatrixP(example.Pfield.size(), 1);

        PtrList<volVectorField> ULmodes;

        for (label k = 0; k < example.liftfield.size(); k++)
        {
            ULmodes.append(example.liftfield[k]);
        }   

        for (label k = 0; k < NmodesUproj; k++)
        {
            ULmodes.append(example.Umodes[k]);
        }


        Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Ufield,
                                 ULmodes, (NmodesUproj + example.liftfield.size()));
        Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Pfield, example.Pmodes,
                                 NmodesPproj );

        PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
               ULmodes, coeffU,(NmodesUproj+ example.liftfield.size()));
        PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                example.Pmodes, coeffP, NmodesPproj );

        Eigen::MatrixXd L2errorProjU = ITHACAutilities::error_listfields(example.Ufield,
                                       rec_fieldU);
        Eigen::MatrixXd L2errorProjP = ITHACAutilities::error_listfields(example.Pfield,
                                       rec_fieldP);

	
        L2errorProjMatrixU.col(0) = L2errorProjU;
        L2errorProjMatrixP.col(0) = L2errorProjP;

	

    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixP, "L2errorProjMatrixP", "eigen",
                               "./ITHACAoutput/l2error");

 
// Resize the modes for projection
    example.Pmodes.resize(NmodesPproj);
    example.Umodes.resize(NmodesUproj);
  
 
if (example.ExplicitMethod == "A"|| example.ExplicitMethod == "B")
{
		Eigen::MatrixXd L2errorProjMatrixPhi(example.Ufield.size(), 1);
	Eigen::MatrixXd coeffPhi = ITHACAutilities::get_coeffs(example.Phifield, example.Phimodes,
                                 NmodesUproj );
	PtrList<surfaceScalarField> rec_fieldPhi = ITHACAutilities::reconstruct_from_coeff(
                example.Phimodes, coeffPhi, NmodesUproj );
	Eigen::MatrixXd L2errorProjPhi = ITHACAutilities::error_listfields(example.Phifield,
                                       rec_fieldPhi);
	L2errorProjMatrixPhi.col(0) = L2errorProjPhi;
    ITHACAstream::exportMatrix(L2errorProjMatrixPhi, "L2errorProjMatrixPhi", "eigen",
                               "./ITHACAoutput/l2error");
  example.Phimodes.resize(NmodesUproj);
}


    // Galerkin Projection
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);

    reducedUnsteadyNSExplicit reduced(example);
    // Set values of the reduced stuff
    reduced.nu = example.mu(0, 0);
    reduced.tstart = example.startTime;
    reduced.finalTime = example.finalTime;
    reduced.dt = example.timeStep;
    reduced.storeEvery = example.timeStep;
    reduced.exportEvery = example.timeStep;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1.0;
    auto start_ROM = std::chrono::high_resolution_clock::now(); 
    reduced.solveOnline(vel_now, 1);
    auto finish_ROM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
    // Reconstruct the solution and export it
    reduced.reconstruct("./ITHACAoutput/ReconstructionExplicit/");
    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields(
                                         example.Pfield, reduced.PREC);
 
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                                "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                                "./ITHACAoutput/l2error");




    exit(0);
}


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
    YJunction.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "reducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorialY: public unsteadyNS
{
    public:
        explicit tutorialY(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
		    U = U0;
		    p = P0;
                    mu_now[0] = mu(0, 0);
                    truthSolve(mu_now);
	            restart();
            }
        }

	void onlineSolveFull( fileName folder)
        {
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (ITHACAutilities::check_folder(folder))
            {
		ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Pfield_on, p, folder);
            }
            else
            {
                mkDir(folder);
                ITHACAutilities::createSymLink(folder);
                label k = para_set_BC;

		U = U0;
		p = P0;
                mu_now[0] = mu(0, 0);
                truthSolve(mu_now, folder);
	        restart();
            }
        }

// Method to compute the lifting functions for this tutorial
void liftSolve()
{
    for (label k = 0; k < 1; k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[0] = 1; // velocity magnitude = 1  and inlet 1
        v1[1] = -1; // velocity magnitude = 1  and inlet 1
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
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

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



    for (label k = 1; k < 2; k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[0] = 1; // velocity magnitude = 1  and inlet 2
        v1[1] = 1; // velocity magnitude = 1  and inlet 2
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
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

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
    // Construct the tutorialY object
    tutorialY example(argc, argv);

// the offline samples for the boundary conditions
    word par_offline_BC("./timeBCoffline");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    Eigen::MatrixXd par_on_BC;
    if (example.bcMethod == "lift")
    {
	word par_online_BC("./timeBConControl");
    	par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    }
    else if (example.bcMethod == "penalty")
    {
	word par_online_BC("./timeBConPenalty");
    	par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    }

    // Read parameters from ITHACAdict file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    int NmodesOut = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 20);
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

    if (example.bcMethod == "lift")
    {
        example.inletIndex.resize(2, 2); 
        example.inletIndex(0, 0) = 1;  // Patch inlet 1
        example.inletIndex(0, 1) = 0;  
        example.inletIndex(1, 0) = 2;  // Patch inlet 2
        example.inletIndex(1, 1) = 0;  
   }
   else if (example.bcMethod == "penalty")
   {
        // Set the inlet boundaries where we have non homogeneous boundary conditions
        example.inletIndex.resize(4, 2); // rows: total number of patches 
        example.inletIndex(0, 0) = 1;  // Patch inlet 1
        example.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
        example.inletIndex(1, 0) = 1;  // Patch inlet 1: y-direction
        example.inletIndex(1, 1) = 1;  // Patch inlet 2
        example.inletIndex(2, 0) = 2;  // Patch inlet 2: x-direction
        example.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
        example.inletIndex(3, 0) = 2;  // Patch inlet 2: x-direction
        example.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
    }

    // Time parameters
    example.startTime = 0;
    example.finalTime = 12;
    example.timeStep = 0.0005;
    example.writeEvery = 0.03;
    example.Dim = 2;

    // Perform The Offline Solve;
    example.offlineSolve();

    // Perform POD
    auto start_POD = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "lift")
    {
            // Search the lift function
   	    example.liftSolve();
  	    // Normalize the lifting function
	    ITHACAutilities::normalizeFields(example.liftfield);
    	    // Create homogeneous basis functions for velocity
    	    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    	    // Perform a POD decomposition for velocity and pressure
   	    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    	    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    }
    else if (example.bcMethod == "penalty")
    {
            // Perform a POD decomposition for velocity and pressure
            ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
            ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    }

    auto finish_POD = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
    std::cout << "elapsed_POD: " << elapsed_POD.count() << " seconds.";
    std::cout << std::endl;

    // PPE method
    NmodesSUPproj = 0;

    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl; 
  
    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.01;
    reduced.tstart = 0;
    reduced.finalTime = 18;
    reduced.dt = 0.0005;

    Eigen::MatrixXd vel_now = par_on_BC;

    reduced.maxIter = 100;
    reduced.tolerance = 1e-5;
    reduced.timeSteps = 5;

    auto start_penalty = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "penalty")
    {
            // set initial quess for penalty factors
	    reduced.tauInit = Eigen::MatrixXd::Zero(4,1);
            reduced.tauInit <<  1e-6,1e-6,1e-6,1e-6; 
            reduced.tauU = reduced.penalty_PPE_time(vel_now, reduced.tauInit);
	
    }
    auto finish_penalty = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_penalty = finish_penalty - start_penalty;
    std::cout << "elapsed_penalty: " << elapsed_penalty.count() << " seconds.";
    std::cout << std::endl;  

    // Set the online temperature BC and solve reduced model
    auto start_ROM = std::chrono::high_resolution_clock::now();
    reduced.solveOnline_PPE(vel_now);
    auto finish_ROM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
    std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
    std::cout << std::endl;

    auto start_ROM_REC = std::chrono::high_resolution_clock::now();
    reduced.reconstruct_sup("./ITHACAoutput/ReconstructionPPE", 60);
    auto finish_ROM_REC = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM_REC = finish_ROM_REC - start_ROM_REC;
    std::cout << "elapsed_ROM_REC: " << elapsed_ROM_REC.count() << " seconds.";
    std::cout << std::endl;
   
    // Performing full order simulation for second parameter set 
    tutorialY HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.01;
    HFonline2.mu_range(0, 1) = 0.01;
    HFonline2.genEquiPar();
    HFonline2.inletIndex.resize(4, 2); // rows: total number of patches 
    HFonline2.inletIndex(0, 0) = 1;  // Patch inlet 1
    HFonline2.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
    HFonline2.inletIndex(1, 0) = 1;  // Patch inlet 1: y-direction
    HFonline2.inletIndex(1, 1) = 1;  // Patch inlet 2
    HFonline2.inletIndex(2, 0) = 2;  // Patch inlet 2: x-direction
    HFonline2.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
    HFonline2.inletIndex(3, 0) = 2;  // Patch inlet 2: x-direction
    HFonline2.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 18;
    HFonline2.timeStep = 0.0005;
    HFonline2.writeEvery = 0.03;
    HFonline2.Dim = 2;

    // Set values BCs for each timeSetp:
    HFonline2.timeBCoff = vel_now;

    // Reconstruct the online solution
    HFonline2.onlineSolveFull(par_on_BC, 1,
                             "./ITHACAoutput/HFonline2");

    // Reading high-fidelity solutions for the parameter set
    // for which the offline solve has been performed (skipping IC)

    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields(
                                         example.Pfield_on, reduced.PREC);
    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2error");


    // Post-Process
    Eigen::MatrixXd PostP(example.Ufield_on.size(), 2); 
    for (label i = 0; i < example.Ufield_on.size(); i++)
    {
	PostP(i, 0) =  0.5*fvc::domainIntegrate(example.Ufield_on[i] & example.Ufield_on[i]).value();

	PostP(i, 1) = 0.5*fvc::domainIntegrate(reduced.UREC[i] & reduced.UREC[i]).value();
	
    }
    ITHACAstream::exportMatrix(PostP, "PostP", "eigen", "./ITHACAoutput/PostProcess");

exit(0); 
} 



  



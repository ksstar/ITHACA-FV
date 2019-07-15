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
    elbowLam.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "reducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorialELB: public unsteadyNS
{
    public:
        explicit tutorialELB(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            Vector<double> inl(0, 0, 0);
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
                for (label i = 0; i < par_BC.cols(); i++)
                {
		    U = U0;
		    p = P0;
                    inl[0] = par_BC(0, i);
                    mu_now[0] = mu(0, 0);
                    assignBC(U, inletIndex(0,0), inl);
                    truthSolve(mu_now);
	            restart();
                }
            }
        }



	void onlineSolveFull(Eigen::MatrixXd par_BC, label para_set_BC, fileName folder)
        {
 	    Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (ITHACAutilities::check_folder(folder))
            {
            }
            else
            {
                mkDir(folder);
                ITHACAutilities::createSymLink(folder);
                label i = para_set_BC;

	        for (label i = 0; i < par_BC.cols(); i++)
                {
		    U = U0;
		    p = P0;
                    inl[0] = par_BC(0, i);
                    mu_now[0] = mu(0, 0);
                    assignBC(U, inletIndex(0,0), inl);
                    truthSolve(mu_now, folder);
	            restart();
                }

            }
        }

        void onlineSolveRead(fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Pfield_on, p, folder);
            }
            else
            {
            }
        }


// Method to compute the lifting function
void liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
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
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        while (simple.correctNonOrthogonal())
        {
                    fvVectorMatrix UEqn
                    (
                        fvm::div(phi, Ulift)
                    );
                    UEqn.solve();
       }

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
    // Construct the tutorialLID object
    tutorialELB example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline("./par_offline");
    Eigen::MatrixXd par_off = ITHACAstream::readMatrix(par_offline);
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online("./par_online");
    Eigen::MatrixXd par_on = ITHACAstream::readMatrix(par_online);
    // Read parameters from ITHACAdict file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 20);
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
    example.inletIndex.resize(2, 2);
    example.inletIndex(0, 0) = 1;
    example.inletIndex(0, 1) = 0;
    example.inletIndex(1, 0) = 2;
    example.inletIndex(1, 1) = 1;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 20;
    example.timeStep = 0.05;
    example.writeEvery = 0.05;

    // Perform The Offline Solve;
    example.offlineSolve(par_off);

    auto start_POD = std::chrono::high_resolution_clock::now();

    if (example.bcMethod == "lift")
    {
	// Solve the supremizer problem
    	//example.solvesupremizer();
        // Search the lift function
   	example.liftSolve();
    	// Create homogeneous basis functions for velocity
    	example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    	// Perform a POD decomposition for velocity and pressure
   	ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    	ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
       // ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
                  //      example.supex, 1, NmodesSUPout);


    }
    else if (example.bcMethod == "penalty")
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesOut);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesOut);
    }

    auto finish_POD = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
    std::cout << "elapsed_POD: " << elapsed_POD.count() << " seconds.";
    std::cout << std::endl;

    // Solve the supremizer problem
    auto start_sup = std::chrono::high_resolution_clock::now();
    example.solvesupremizer("modes");
    auto finish_sup = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sup = finish_sup - start_sup;
    std::cout << "elapsed_sup: " << elapsed_sup.count() << " seconds.";
    std::cout << std::endl;

    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl;

    // Resize the modes for projection
    example.Umodes.resize(NmodesUproj);
    example.Pmodes.resize(NmodesPproj);

    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.01;
    reduced.tstart = 0;
    reduced.finalTime = 20;
    reduced.dt = 0.05;
    reduced.maxIter = 100;
    reduced.tolerance = 0.0001;
    reduced.timeSteps = 10;

    auto start_ROM = std::chrono::high_resolution_clock::now();
    // Set the online temperature BC and solve reduced model
    //for (label k = 0; k < (par_on.rows()); k++)
    //{
        Eigen::MatrixXd vel_now(2, 1);
        vel_now(0, 0) = 1;
        vel_now(1, 0) = 3;
        //vel_now(1, 0) = par_on(k, 1);

        if (example.bcMethod == "penalty")
        {
            // set initial quess for penalty factors temperature BCs (>  0)
            Eigen::MatrixXd tauInit(1,1) ;
	    ///Eigen::MatrixXd tauInit(vel_now.rows(),1) ;
            tauInit << 1.0;
            reduced.tauU = reduced.penalty_sup(vel_now, tauInit);
	    reduced.tauU = tauInit;
        }

        reduced.solveOnline_sup(vel_now, 0, 2);
	//reduced.solveOnline_sup(vel_now, k, par_on.rows());
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 1);
    //}
    auto finish_ROM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
    std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
    std::cout << std::endl;

exit(0);
}





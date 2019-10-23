/*---------------------------------------------------------------------------*\
     ???????????????  ??? ??????  ??????? ??????       ???????????   ???
     ???????????????  ???????????????????????????      ???????????   ???
     ???   ???   ???????????????????     ????????????????????  ???   ???
     ???   ???   ???????????????????     ????????????????????  ???? ????
     ???   ???   ???  ??????  ??????????????  ???      ???      ???????
     ???   ???   ???  ??????  ??? ??????????  ???      ???       ?????

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
    Example of an unsteady NS Reduction Problem with time-dependent boundary
    conditions.
SourceFiles
    17YJunction.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorial17: public unsteadyNS
{
    public:
        explicit tutorial17(int argc, char* argv[])
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
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/",0,2);
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/",0,2);
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

};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial17 object
    tutorial17 example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./timeBCoff");
    example.timeBCoff = ITHACAstream::readMatrix(par_offline_BC);
      Eigen::MatrixXd par_on_BC;
    word par_online_BC("./timeBCon");
    par_on_BC = ITHACAstream::readMatrix(par_online_BC);

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
    example.mu_range(0, 0) = 0.0000002;
    example.mu_range(0, 1) = 0.0000002;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();

    // Set the inlet boundaries where we have non homogeneous VELOCITY boundary conditions
    example.inletIndex.resize(1, 2); // rows: total number of patches
    example.inletIndex(0, 0) = 2;  // Patch 2 - Inlet
    example.inletIndex(0, 1) = 2;  // Patch 2 - Inlet : z-direction
    example.inletPatch.resize(1, 1);
    example.inletPatch(0, 0) = example.inletIndex(0, 0);  // Patch 2 Inlet

    // Set the inlet boundaries where we have non homogeneous PRESSURE boundary conditions
    example.inletIndexP.resize(1, 2); // rows: total number of patches
    example.inletIndexP(0, 0) = 1;  // Patch 1 - Inlet

    // Time parameters
    example.startTime = 0;
    example.finalTime = 2000;
    example.timeStep = 5;
    example.writeEvery = 5;

    // Perform The Offline Solve;
    example.offlineSolve();

    // Perform POD
    auto start_POD = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "penalty" || example.bcMethod == "penaltyP")
    {
        // Perform a POD decomposition for velocity and pressure
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                     NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                    NmodesPout);
	//ITHACAPOD::getModesVelBasis(example.Pfield, example.Pmodes, example.Ufield,  example.podex, 0, 0,
          //           NmodesPout);
	// Solve the supremizer problem
  	 //example.solvesupremizer("modes");
    }

    auto finish_POD = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
    std::cout << "elapsed_POD: " << elapsed_POD.count() << " seconds.";
    std::cout << std::endl;

    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl;

    reducedUnsteadyNS reduced(example);
    // Set values of the online phase
    reduced.nu = 0.0000002;;
    reduced.tstart = 0;
    reduced.finalTime = 2000;
    reduced.dt = 5;
    reduced.storeEvery = 5;
    reduced.exportEvery = 5;
    // Set values velocity boundary conditions of the online phase
    Eigen::MatrixXd vel_now = par_on_BC;
    // Set values of the iterative penalty method
    reduced.maxIterPenalty = 100;
    reduced.tolerancePenalty = 1e-5;
    reduced.timeStepPenalty = 5;
 
    auto start_penalty = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "penaltyP")
    {
            // Set initial quess for penalty factors
	    reduced.tauIter = Eigen::MatrixXd::Zero(1,1);
            reduced.tauIter <<  0.01;
	    // Solve for the penalty factors with the iterative solver
            reduced.taup = reduced.tauIter;
//reduced.penalty_PPE(vel_now, reduced.tauIter);
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
    reduced.reconstruct_PPE("./ITHACAoutput/ReconstructionPPE");
    auto finish_ROM_REC = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM_REC = finish_ROM_REC - start_ROM_REC;
    std::cout << "elapsed_ROM_REC: " << elapsed_ROM_REC.count() << " seconds.";
    std::cout << std::endl;

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::error_listfields(
                                         example.Pfield, reduced.PREC);

//Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");

exit(0);
}







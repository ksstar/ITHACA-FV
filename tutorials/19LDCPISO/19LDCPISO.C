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
    19LidDrivenCavityClean.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNSPISO.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNSPISO.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial19PISO: public unsteadyNSPISO
{
    public:
        explicit tutorial19PISO(int argc, char* argv[])
            :
            unsteadyNSPISO(argc, argv),
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
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial19PISO object
    tutorial19PISO example(argc, argv);
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
    example.finalTime = 1.0;
    example.timeStep = 0.01;
    example.writeEvery = 0.01;
    // Perform The Offline Solve;
    auto start_FOM = std::chrono::high_resolution_clock::now();
    example.offlineSolve();
    auto finish_FOM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_FOM = finish_FOM - start_FOM;
    // Solve the supremizer problem
    example.solvesupremizer();
    if (example.bcMethod == "lift")
    {
    // Search the lift function
    example.liftSolve();
    // Normalize the lifting function
    //ITHACAutilities::normalizeFields(example.liftfield);
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    }
  
    
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
                        example.supex, 1, NmodesSUPout);
/*
     Eigen::MatrixXd coeffL2Uortho = ITHACAutilities::get_coeffs_ortho(example.Ufield,
              example.Umodes, NmodesUout);
    ITHACAstream::exportMatrix(coeffL2Uortho, "coeffL2Uortho", "eigen",
                               "./ITHACAoutput/Matrices/");

Eigen::MatrixXd coeffL2U = ITHACAutilities::get_coeffs(example.Ufield,
              example.Umodes, NmodesUout);
    ITHACAstream::exportMatrix(coeffL2U, "coeffL2U", "eigen",
                               "./ITHACAoutput/Matrices/");

Eigen::MatrixXd coeffL2Portho = ITHACAutilities::get_coeffs_ortho(example.Pfield,
              example.Pmodes, NmodesPout);
    ITHACAstream::exportMatrix(coeffL2Portho, "coeffL2Portho", "eigen",
                               "./ITHACAoutput/Matrices/");

Eigen::MatrixXd coeffL2P = ITHACAutilities::get_coeffs_ortho(example.Pfield,
              example.Pmodes, NmodesPout);
    ITHACAstream::exportMatrix(coeffL2P, "coeffL2P", "eigen",
                               "./ITHACAoutput/Matrices/");*/


    // Galerkin Projection
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);

 /* //Post-Process
    Eigen::MatrixXd PostPfom(example.Ufield.size(), 1); 

    for (label i = 0; i < example.Ufield.size(); i++)
    {
	PostPfom(i, 0) =  0.5*fvc::domainIntegrate(example.Ufield[i] & example.Ufield[i]).value();
    }

    ITHACAstream::exportMatrix(PostPfom, "PostP", "eigen", "./ITHACAoutput/PostProcess");*/
 
   reducedUnsteadyNSPISO reduced(example);
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
    reduced.solveOnline_sup(vel_now, 1);
    auto finish_ROM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
    // Reconstruct the solution and export it
    reduced.reconstruct_sup("./ITHACAoutput/ReconstructionPISO/");
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

 /*   //Post-Process
    Eigen::MatrixXd PostP(example.Ufield.size(), 2); 

    for (label i = 0; i < example.Ufield.size(); i++)
    {
	PostP(i, 0) =  0.5*fvc::domainIntegrate(example.Ufield[i] & example.Ufield[i]).value();

	PostP(i, 1) = 0.5*fvc::domainIntegrate(reduced.UREC[i] & reduced.UREC[i]).value();
	
    }

    ITHACAstream::exportMatrix(PostP, "PostP", "eigen", "./ITHACAoutput/PostProcess");*/

    //Computational Time
    std::cout << "elapsed_FOM: " << elapsed_FOM.count() << " seconds.";
    std::cout << std::endl;
    std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
    std::cout << std::endl;





    exit(0);
}


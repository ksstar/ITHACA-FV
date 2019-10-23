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
    LidDrivenCavity.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorialLID2: public unsteadyNS
{
    public:
        explicit tutorialLID2(int argc, char* argv[])
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
                for (label i = 0; i < par_BC.rows(); i++)
                {
		    U = U0;
		    p = P0;
                    inl[inletIndex(0,i+1)] = par_BC(i, 0);
                    mu_now[0] = mu(0, 0);
                    assignBC(U, inletIndex(i,0), inl);
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
                label k = para_set_BC;

	        for (label i = 0; i < par_BC.rows(); i++)
                {
		    U = U0;
		    p = P0;
                    inl[inletIndex(0,i+1)] = par_BC(i, k);
                    mu_now[0] = mu(0, 0);
                    assignBC(U, inletIndex(i,0), inl);
                    truthSolve(mu_now, folder);
	            restart();
                }

            }
        }

};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorialLID object
    tutorialLID2 example(argc, argv);
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
    example.mu_range(0, 0) = 0.0001;
    example.mu_range(0, 1) = 0.0001;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 10;
    example.timeStep = 0.0005;
    example.writeEvery = 0.01;

    // Perform The Offline Solve;
    example.offlineSolve(par_off);

    auto start_POD = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "lift")
    {

            // Search the lift function
   	    example.liftSolve();
	    // Normalize the lifting function
	   // ITHACAutilities::normalizeFields(example.liftfield);
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

//exit(0);

 /*  // Solve the supremizer problem
    auto start_sup = std::chrono::high_resolution_clock::now();
    example.solvesupremizer("modes");
    auto finish_sup = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sup = finish_sup - start_sup;
    std::cout << "elapsed_sup: " << elapsed_sup.count() << " seconds.";
    std::cout << std::endl; */


    // PPE method
    NmodesSUPproj = 0;

    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl; 


/* // Create a list with number of modes for which the projection needs to be performed
    Eigen::MatrixXd List_of_modes(NmodesOut, 1);
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        List_of_modes(i, 0) = i + 1;
    }

    // Export with number of modes for which the projection needs to be performed
    ITHACAstream::exportMatrix(List_of_modes, "List_of_modes", "eigen", "./ITHACAoutput/l2error");

    PtrList<volVectorField> Umodes;

    if (example.liftfield.size() != 0)
    {
        for (label k = 0; k < example.liftfield.size(); k++)
        {
            Umodes.append(example.liftfield[k]);
        }
    }

    if (example.aveMethod == "mean")
    {
      Umodes.append(example.Umean[0]);
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
            Umodes.append(example.Umodes[k]);
    }

    if (example.supex == 1)
    {
        for (label k = 0; k < NmodesSUPproj; k++)
        {
            Umodes.append(example.supmodes[k]);
        }
    }

    PtrList<volScalarField> Pmodes;

    if (example.aveMethod == "mean")
    {
       Pmodes.append(example.Pmean[0]);
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        Pmodes.append(example.Pmodes[k]);
    }
   


    // Perform the projection for all number of modes in list List_of_modes
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size()-1, List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixP(example.Pfield.size()-1, List_of_modes.rows());

    // Calculate the coefficients and L2 error and store the error in a matrix for each number of modes
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
	std::cout << "Number of modes = " << i << std::endl;

	if (example.aveMethod == "mean")
        { 
	    Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Usub,Umodes,
                                   List_of_modes(i, 0) + example.liftfield.size() + example.Umean.size() + NmodesSUPproj);
            Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Psub, Pmodes,
                                 List_of_modes(i, 0) + example.Pmean.size());
            PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                		 Umodes, coeffU, List_of_modes(i, 0) + example.liftfield.size() + example.Umean.size() + NmodesSUPproj);
            PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                		 Pmodes, coeffP, List_of_modes(i, 0) + example.Pmean.size());

            L2errorProjMatrixU.col(i) = ITHACAutilities::error_listfields_min_IC(example.Usub, rec_fieldU);
            L2errorProjMatrixP.col(i) = ITHACAutilities::error_listfields_min_IC(example.Psub, rec_fieldP);

	}
	else
	{
            Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Ufield,Umodes,
                                   List_of_modes(i, 0) + example.liftfield.size() + NmodesSUPproj);
            Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Pfield, Pmodes,
                                 List_of_modes(i, 0));
            PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                		 Umodes, coeffU, List_of_modes(i, 0));
            PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                		 Pmodes, coeffP, List_of_modes(i, 0));
            L2errorProjMatrixU.col(i) = ITHACAutilities::error_listfields_min_IC(example.Ufield, rec_fieldU);
            L2errorProjMatrixP.col(i) = ITHACAutilities::error_listfields_min_IC(example.Pfield, rec_fieldP);
	}
    } 

    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixP, "L2errorProjMatrixP", "eigen",
                               "./ITHACAoutput/l2error"); 

    // Resize the modes for projection
    example.Umodes.resize(NmodesUproj);
    example.Pmodes.resize(NmodesPproj); */

    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.0001;
    reduced.tstart = 0;
    reduced.finalTime = 10;
    reduced.dt = 0.0005;
    // Set values of the iterative penalty method
    reduced.maxIterPenalty = 100;
    reduced.tolerancePenalty = 1e-5;
    reduced.timeStepPenalty = 5;
    reduced.storeEvery = 0.01;
    reduced.exportEvery = 0.01;

    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;

    auto start_penalty = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "penalty")
    {
            // set initial quess for penalty factors
	    reduced.tauIter = Eigen::MatrixXd::Zero(1,1);
            reduced.tauIter << 0.01;
            //reduced.tauU = reduced.penalty_PPE(vel_now, reduced.tauInit);
            reduced.tauU = reduced.tauIter;
    }

    auto finish_penalty = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_penalty = finish_penalty - start_penalty;
    std::cout << "elapsed_penalty: " << elapsed_penalty.count() << " seconds.";
    std::cout << std::endl;  

//exit(0);
    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (1); k++)
    {
	auto start_ROM = std::chrono::high_resolution_clock::now();
        vel_now(0, 0) = 1;        
        reduced.solveOnline_PPE(vel_now, k);
        auto finish_ROM = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
        std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
        std::cout << std::endl;

//exit(0);
        auto start_ROM_REC = std::chrono::high_resolution_clock::now();
        reduced.reconstruct_PPE("./ITHACAoutput/ReconstructionPPE");
        auto finish_ROM_REC = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM_REC = finish_ROM_REC - start_ROM_REC;
        std::cout << "elapsed_ROM_REC: " << elapsed_ROM_REC.count() << " seconds.";
        std::cout << std::endl;
    }

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

//Post-Process
    Eigen::MatrixXd PostP(example.Ufield.size(), 2); 

    for (label i = 0; i < example.Ufield.size(); i++)
    {
	PostP(i, 0) = 0.5*fvc::domainIntegrate(example.Ufield[i] & example.Ufield[i]).value();
	PostP(i, 1) = 0.5*fvc::domainIntegrate(reduced.UREC[i] & reduced.UREC[i]).value();	
    }

    ITHACAstream::exportMatrix(PostP, "PostP", "eigen", "./ITHACAoutput/PostProcess");


exit(0); 
} 



  



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
    Example of steady NS Reduction Problem
SourceFiles
    03steadyNS.C
\*---------------------------------------------------------------------------*/

#include "steadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "reducedSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"

class tutorial03_UQ : public steadyNS
{
    public:
        /// Constructor
        explicit tutorial03_UQ(int argc, char* argv[])
            :
            steadyNS(argc, argv),
            U(_U()),
            p(_p()),
            args(_args())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        /// Arg List
        argList& args;

        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else
            {
                Vector<double> Uinl(0, 0, 0);

                for (label i = 0; i < mu.cols(); i+=4)
                {
                    U = U0;
                    p = P0;
                    mu_now[0] = mu(0, i);
                    change_viscosity(mu(0, i));
                    truthSolve(mu_now);
                    restart();
                }
            }
        }


	/// Perform an Offline solve
        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                //mu_samples =
                  //  ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else
            {
                Vector<double> Uinl(0, 0, 0);
                for (label j = 0; j < par_BC.cols(); j++)
                {
                    U = U0;
                    p = P0;
                    mu_now[0] = mu(0, j);
                    change_viscosity(mu(0, j));
	            Uinl[0] = par_BC(0, j);
                    //Uinl[0] = 1;
                    assignBC(U, inletIndex(0, 0), Uinl);
                    truthSolve(mu_now);
                    restart();
                }
            }
        }
	/// Perform an Offline solve
        void offlineSolve2(Eigen::MatrixXd par_BC, word folder)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;

	    if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline2/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline2/");
                //mu_samples =
                  //  ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
	    else
            {
            
                Vector<double> Uinl(0, 0, 0);
                for (label j = 0; j < par_BC.cols(); j++)
                {
                    U = U0;
                    p = P0;
                    mu_now[0] = mu(0, j);
                    change_viscosity(mu(0, j));
	            Uinl[0] = par_BC(0, j);
                    //Uinl[0] = 1;
                    assignBC(U, inletIndex(0, 0), Uinl);
                    truthSolve(mu_now, folder );
                    restart();
                }
	     }
            
        }
};


int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial03_UQ example(argc, argv);
    // the offline samples 
    Info << "Here" << endl;
    word par_offline_BC("./par_offline_vel");
    word par_offline_nu("./par_offline_nu");
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online_BC("./par_online_vel");
    word par_online_nu("./par_online_nu");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    Eigen::MatrixXd par_off_nu = ITHACAstream::readMatrix(par_offline_nu);
    Eigen::MatrixXd par_on_nu = ITHACAstream::readMatrix(par_online_nu);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    // Read the par file where the parameters are stored
    //word filename("./par_online");
    example.mu = ITHACAstream::readMatrix(par_online_nu);
    // Set the number of parameters
    //example.Pnumber = 1;
    // Set samples
    //example.Tnumber = 1;
    // Set the parameters infos
    //example.setParameters();
    // Set the parameter ranges
    //example.mu_range(0, 0) = 1;
    //example.mu_range(0, 1) = 1;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    //example.offlineSolve();
    example.offlineSolve(par_off_BC);	

    auto start_sup = std::chrono::high_resolution_clock::now();
    // Solve the supremizer problem
    example.solvesupremizer();
    auto finish_sup = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sup = finish_sup - start_sup;

    auto start_POD = std::chrono::high_resolution_clock::now();
    // Search the lift function for the velocity
    example.liftSolve();
    // Create homogeneous basis functions for temperature
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    //ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity pressure and supremizers and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
                        example.supex, 1, NmodesSUPout);
    auto finish_POD = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
// Create a list with number of modes for which the projection needs to be performed
    Eigen::MatrixXd List_of_modes(NmodesOut, 1);

    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        List_of_modes(i, 0) = i + 1;
    }

    // Export with number of modes for which the projection needs to be performed
    ITHACAstream::exportMatrix(List_of_modes, "List_of_modes", "eigen",
                               "./ITHACAoutput/l2error");

    // Create locally the velocity modes
    PtrList<volVectorField> ULmodes;

    for (label k = 0; k < example.liftfield.size(); k++)
    {
        ULmodes.append(example.liftfield[k]);
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        ULmodes.append(example.Umodes[k]);
    }

    for (label k = 0; k < NmodesSUPproj; k++)
    {
        ULmodes.append(example.supmodes[k]);
    }

    // Perform the projection for all number of modes in list List_of_modes
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size(), List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixP(example.Pfield.size(), List_of_modes.rows());

    // Calculate the coefficients and L2 error and store the error in a matrix for each number of modes
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Ufield, example.Umodes,
                                 List_of_modes(i, 0) + example.liftfield.size());
        Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Pfield, example.Pmodes, List_of_modes(i, 0));
        PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                example.Umodes, coeffU, List_of_modes(i, 0));
        PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                example.Pmodes, coeffP, List_of_modes(i, 0));
        Eigen::MatrixXd L2errorProjU = ITHACAutilities::error_listfields(example.Ufield,
                                       rec_fieldU);
        Eigen::MatrixXd L2errorProjP = ITHACAutilities::error_listfields(example.Pfield,
                                       rec_fieldP);
        L2errorProjMatrixU.col(i) = L2errorProjU;
        L2errorProjMatrixP.col(i) = L2errorProjP;
    }

    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixP, "L2errorProjMatrixP", "eigen",
                               "./ITHACAoutput/l2error");

    auto start_matrix = std::chrono::high_resolution_clock::now();
    // Perform the Galerkin Projection
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;

// Resize the modes for projection
    example.Pmodes.resize(NmodesPproj);
    example.Umodes.resize(NmodesUproj);
    // Create the reduced object
    reducedSteadyNS ridotto(example);
    // Set the inlet velocity
    Eigen::MatrixXd vel_now(2, 1);
    //vel_now(0, 0) = 1;
    //vel_now(1, 0) = 0;

    auto start_ROM = std::chrono::high_resolution_clock::now();
    // Perform an online solve for the new values of inlet velocities
    //for (label k = 0; k <  par_on_BC.cols(); k++)
    for (label k = 0; k <  par_on_nu.cols(); k++)
    {
        // Set the reduced viscosity
        ridotto.nu = example.mu(0, k);
	vel_now(0, 0) = par_on_BC(0,k);
        //vel_now(0, 0) = 1;
        ridotto.solveOnline_sup(vel_now);
        Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
        ridotto.online_solution.append(tmp_sol);
    }
    

    
    // Reconstruct and export the solution
    ridotto.reconstruct_sup("./ITHACAoutput/Reconstruction/");
    auto finish_ROM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;

    // Save the online solution
    ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");


    tutorial03_UQ HFonline2(argc, argv);
    HFonline2.mu = ITHACAstream::readMatrix(par_online_nu);
    HFonline2.inletIndex.resize(1, 2);
    HFonline2.inletIndex(0, 0) = 0;
    HFonline2.inletIndex(0, 1) = 0;

    auto start_FOM = std::chrono::high_resolution_clock::now();
    // Perform the offline solve
    HFonline2.offlineSolve2(par_on_BC, "./ITHACAoutput/Offline2/");
    auto finish_FOM = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_FOM = finish_FOM - start_FOM;

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         HFonline2.Ufield, ridotto.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields(
                                         HFonline2.Pfield, ridotto.PREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2error");

    //Post-Process
    Eigen::MatrixXd PostP(HFonline2.Ufield.size(), 6); 

    for (label i = 0; i < HFonline2.Ufield.size(); i++)
    {
	PostP(i, 0) =  0.5*fvc::domainIntegrate(HFonline2.Ufield[i] & HFonline2.Ufield[i]).value();
	
	// min/max Outlet (patch 1) velocity
	PostP(i, 1) = max(HFonline2.Ufield[i].boundaryField()[1].component(0) );
        PostP(i, 2) = min(HFonline2.Ufield[i].boundaryField()[1].component(0) );

	PostP(i, 3) = 0.5*fvc::domainIntegrate(ridotto.UREC[i] & ridotto.UREC[i]).value();
	PostP(i, 4) = max(ridotto.UREC[i].boundaryField()[1].component(0) );
        PostP(i, 5) = min(ridotto.UREC[i].boundaryField()[1].component(0) );
	

    }

    ITHACAstream::exportMatrix(PostP, "PostP", "eigen",
                               "./ITHACAoutput/PostProcess");



        
                
          
     




    // Post-processing


    return 0;
}

//--------
/// \dir 03steadyNS Folder of the turorial 3
/// \file
/// \brief Implementation of a tutorial of a steady Navier-Stokes problem

/// \example 03steadyNS.C
/// \section intro_sreadyNS Introduction to tutorial 3
/// The problems consists of steady Navier-Stokes problem with parametrized viscosity.
/// The physical problem is the backward facing step depicted in the following image:
/// \image html step.png
/// At the inlet a uniform and constant velocity equal to 1 m/s is prescribed.
///
/// \section code03 A detailed look into the code
///
/// In this section are explained the main steps necessary to construct the tutorial N°3
///
/// \subsection header The necessary header files
///
/// First of all let's have a look to the header files that needs to be included and what they are responsible for:
///
/// The header file of ITHACA-FV necessary for this tutorial
///
/// \dontinclude 03steadyNS.C
/// \skip steadyNS
/// \until reducedSteady
///
/// \subsection classtuto03 Implementation of the tutorial03 class
///
/// Then we can define the tutorial03 class as a child of the steadyNS class
/// \skipline tutorial03
/// \until {}
///
/// The members of the class are the fields that needs to be manipulated during the
/// resolution of the problem
///
/// Inside the class it is defined the offlineSolve method according to the
/// specific parametrized problem that needs to be solved.
///
/// \skipline offlineSolve
/// \until {
///
///
/// If the offline solve has already been performed than read the existing snapshots
///
/// \skipline offline
/// \until }
///
/// else perform the offline solve where a loop over all the parameters is performed:
///
/// \skipline else
/// \until }
/// \skipline }
///
/// See also the steadyNS class for the definition of the methods.
///
/// \subsection main Definition of the main function
///
/// Once the tutorial03 class is defined the main function is defined,
/// an example of type tutorial03 is constructed:
///
/// \skipline tutorial03
///
/// In this case the vector of parameter is read from a txt file
///
/// \skipline word
/// \until example.mu
///
/// The inlet boundary is set:
///
/// \skipline example.inlet
/// \until example.inletIndex(0, 1) = 0;
///
/// and the offline stage is performed:
///
/// \skipline Solve()
///
/// and the supremizer problem is solved:
///
/// \skipline supremizer()
///
/// In order to show the functionality of reading fields in this case the lifting function is read
/// from a precomputed simulation with a unitary inlet velocity:
///
/// \skipline stream
///
/// Then the snapshots matrix is homogenized:
///
/// \skipline computeLift
///
/// and the modes for velocity, pressure and supremizers are obtained:
///
/// \skipline getModes
/// \until supfield
///
/// then the projection onto the POD modes is performed with:
///
/// \skipline projectSUP
///
/// the reduced object is constructed:
///
/// \skipline reducedSteady
///
/// and the online solve is performed for some values of the viscosity:
///
/// \skipline Eigen::
/// \until }
///
/// The vel_now matrix in this case is not used since there are no parametrized boundary conditions.
///
/// The viscosity is set with the command:
///
/// \code
/// ridotto.nu = example.mu(k,0)
/// \endcode
///
/// finally the online solution stored during the online solve is exported to file in three different
/// formats with the lines:
///
/// \skipline exportMatrix
/// \until "eigen"
///
/// and the online solution is reconstructed and exported to file
///
/// \skipline reconstruct
///
///
///
///
///
/// \section plaincode The plain program
/// Here there's the plain code
///







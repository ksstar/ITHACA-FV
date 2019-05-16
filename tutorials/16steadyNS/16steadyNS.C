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
    16steadyNS.C
\*---------------------------------------------------------------------------*/

#include "steadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "reducedSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"


class tutorial16 : public steadyNS
{
    public:
        /// Constructor
        explicit tutorial16(int argc, char* argv[])
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

                for (label i = 0; i < mu.cols(); i++)
                {
                    U = U0;
                    p = P0;
                    mu_now[0] = mu(0, i);
                    change_viscosity(mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }

     /// Perform an Offline solve
        void offlineSolve2(word folder)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;

                Vector<double> Uinl(0, 0, 0);

                for (label i = 0; i < mu.cols(); i++)
                {
                    U = U0;
                    p = P0;
                    mu_now[0] = mu(0, i);
                    change_viscosity(mu(0, i));
                    truthSolve(mu_now, folder);
                }
          
        }
};


int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial16 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();
    // Solve the supremizer problem
    example.solvesupremizer();
    ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity pressure and supremizers and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
                        example.supex, 1, NmodesSUPout);
    // Perform the Galerkin Projection
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    // Create the reduced object
    reducedSteadyNS ridotto(example);
    // Set the inlet velocity
    Eigen::MatrixXd vel_now(2, 1);
    vel_now(0, 0) = 1;
    vel_now(1, 0) = 0;

    // Perform an online solve for the new values of inlet velocities
    for (label k = 0; k < example.mu.size(); k++)
    {
        // Set the reduced viscosity
        ridotto.nu = example.mu(0, k);
        ridotto.solveOnline_sup(vel_now);
        Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
        ridotto.online_solution.append(tmp_sol);
    }

    // Save the online solution
    ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");
    // Reconstruct and export the solution
    ridotto.reconstruct_sup("./ITHACAoutput/Reconstruction/");

   // tutorial16 HFonline2(argc, argv);
  //  HFonline2.mu = ITHACAstream::readMatrix(filename);
  //  HFonline2.inletIndex.resize(1, 2);
 //  HFonline2.inletIndex(0, 0) = 0;
  //  HFonline2.inletIndex(0, 1) = 0;

 //   HFonline2.offlineSolve2("./ITHACAoutput/Offline2/");
    
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield, ridotto.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields(
                                         example.Pfield, ridotto.PREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2error");



    return 0;
}




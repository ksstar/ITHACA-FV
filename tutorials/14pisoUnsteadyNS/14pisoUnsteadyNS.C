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

#include "unsteadyNS_piso.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "reducedPisoUnsteadyNS.H"
#include "forces.H"
#include "IOmanip.H"


class tutorial14 : public unsteadyNS_piso
{
    public:
        /// Constructor
        explicit tutorial14(int argc, char* argv[])
            :
            unsteadyNS_piso(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;

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
            }
            else
            {
                Vector<double> Uinl(0.0001, 0, 0);
		scalar pinl = 0;
                label BCind = 0;

                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
			Info << "mu: " << mu(0, i)<< endl;
                    //change_viscosity(mu(0, i));
                    U = U0;
		    p = P0;
                    truthSolve2(mu_now);
                }
            }
        }

};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial14 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    /// Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    /// Time parameters
    example.startTime = 0.0;
    example.finalTime = 0.01;
    example.timeStep = 0.0001;
    example.writeEvery = 0.0001;
    // Perform the Offline Solve;
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();
    // Solve the supremizer problem
    example.solvesupremizer();
    // Search the lift function
    example.liftSolve();
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    // Create the reduced object
    reducedPisoUnsteadyNS reduced(example);
    PtrList<volVectorField> U_rec_list;
    PtrList<volScalarField> P_rec_list;

    // Set values of the online solve
    reduced.nu = 0.01;
    reduced.tstart = 0.0;
    reduced.finalTime = 0.01;
    reduced.dt = 0.0001;
    Eigen::MatrixXd vel(1, 1);
    vel(0, 0) = 1;
    reduced.setOnlineVelocity(vel);
    reduced.solveOnline_Piso();

     //   reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
      //  reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 2);


    exit(0);
}

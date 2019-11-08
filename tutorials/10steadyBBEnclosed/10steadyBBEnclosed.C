/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
    Example of Boussinesq approximation for two way coupling NS-momentum equation
    and heat transport equation for enclosed flows.
SourceFiles
    10UnsteadyBBEnclosed.C
\*---------------------------------------------------------------------------*/

#include "SteadyBB.H"
#include "ITHACAPOD.H"
#include "ReducedSteadyBB.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>
#include "steadyNS.H"
#include "ITHACAPOD.H"
#include "forces.H"

class tutorial10Steady: public SteadyBB
{
    public:
        explicit tutorial10Steady(int argc, char* argv[])
            :
            SteadyBB(argc, argv),
            U(_U()),
            p(_p()),
            p_rgh(_p_rgh()),
            T(_T()),
	    args(_args())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& p_rgh;
        volScalarField& T;

	/// Arg List
        argList& args;

        void offlineSolve(Eigen::MatrixXd par_BC)
        {

            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;
            volScalarField T0 = T;
            volScalarField p0_rgh = p_rgh;


            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
	        ITHACAstream::read_fields(Prghfield, p_rgh, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
            }
            else
            {
		

		//for (label i = 0; i < mu.cols(); i++)
	        for (label i = 0; i < par_BC.rows(); i++)
                {     
		    
		      U = U0;
		      T = T0;
		      p = P0;
		      p_rgh = p0_rgh;

		     for (label j = 0; j < inletIndexT.size(); j++)
                     { 
		         assignBC(T, inletIndexT(j, 0), par_BC(i, j)); 
		     }

                      //mu_now[0] = mu(0, i);
		      //change_viscosity(mu(0, i));

		      truthSolve(mu_now);
		      restart();
                }

            }
        }

        // Method to compute the lifting function for temperature
        void liftSolveT()
        {
            for (label k = 0; k < inletIndexT.rows(); k++)
            {   fvMesh& mesh = _mesh();
                volScalarField T = _T();
                simpleControl simple(mesh);
                volScalarField& alphat = _alphat();
                dimensionedScalar& Pr = _Pr();
                dimensionedScalar& Prt = _Prt();
                label BCind = inletIndexT(k, 0);
                volScalarField Tlift("Tlift" + name(k), T);
                Info << "Solving a lifting Problem" << endl;
                scalar t1 = 1;
                scalar t0 = 0;

                for (label j = 0; j < T.boundaryField().size(); j++)
                {
                    assignIF(Tlift, t0);

                    if (j == BCind)
                    {
                        assignBC(Tlift, j, t1);
                    }
                    else if (T.boundaryField()[BCind].type() == "fixedValue")
                    {
                        assignBC(Tlift, j, t0);
                    }
                    else
                    {
                    }
                }

                while (simple.correctNonOrthogonal())
                {
                    alphat = turbulence->nut() / Prt;
                    alphat.correctBoundaryConditions();
                    volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(alphaEff, Tlift)
                    );
                    TEqn.solve();
                }

                Tlift.write();
                liftfieldT.append(Tlift);
            }
        }
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial10Steady example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./par_offline_BC");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online_BC("./par_online_BC");
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesTproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 100;
    /// Set the parameters infos
    example.setParameters();
    /// Set the parameter ranges
    example.mu_range(0, 0) = 0.00001;
    example.mu_range(0, 1) = 0.00001; //0.0002
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    /// Set the inlet Temperature boundaries where there are non homogeneous boundary conditions
    example.inletIndexT.resize(2, 1);
    example.inletIndexT << 1, 2;
    // Perform the Offline Solve;
    example.offlineSolve(par_off_BC);
    // Search the lift function for the temperature
    example.liftSolveT();
    // Create homogeneous basis functions for temperature
    example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
    // Perform a POD decomposition for velocity temperature and pressure fields
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Prghfield, example.Prghmodes, example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example.podex, 0, 0,
                        NmodesOut);
    // Solve the supremizer problem based on the pressure modes
    example.solvesupremizer("modes");
    // Get reduced matrices
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesTproj,
                       NmodesSUPproj);

    // Create the reduced object
    ReducedSteadyBB reduced(example);
    // Set the inlet velocity
    Eigen::MatrixXd vel(1, 1);
    vel(0, 0) = 0;
    

    // Perform an online solve for the new values of inlet velocities
    //for (int k = 0; k < example.mu.size(); k++)
    for (int k = 0; k < par_on_BC.rows(); k++)
    {
	// Set the temperature boundary conditions
    	Eigen::MatrixXd temp(2, 1);
    	temp(0, 0) = par_on_BC(k, 0);
    	temp(1, 0) = par_on_BC(k, 1);

        // Set the reduced viscosity
        // reduced.nu = example.mu(0, k);
	// std::cout << "nu: " << reduced.nu << std::endl;
	reduced.nu = example.mu(0, 0);
        reduced.Pr = 0.71;
        reduced.solveOnline_sup(temp, vel, k);
        Eigen::MatrixXd tmp_sol(reduced.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(reduced.y.rows()) = reduced.y;
        reduced.online_solution.append(tmp_sol);
    }

    // Save the online solution
    ITHACAstream::exportMatrix(reduced.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");

    reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 1);
    
    return 0;
  
}

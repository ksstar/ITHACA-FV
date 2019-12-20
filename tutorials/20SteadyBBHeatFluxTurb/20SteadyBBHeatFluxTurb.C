/*---------------------------------------------------------------------------*\
Copyright (C) 2019 by the ITHACA-FV authors

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
    20SteadyBBHeatFluxTurb.C
\*---------------------------------------------------------------------------*/

#include "SteadyBBTurb.H"
#include "SteadyBB.H"
#include "steadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedSteadyBBTurb.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>
#include "steadyNS.H"
#include "ITHACAPOD.H"
#include "forces.H"

class tutorial20Steady: public SteadyBBTurb
{
    public:
        explicit tutorial20Steady(int argc, char* argv[])
            :
            SteadyBBTurb(argc, argv),
            U(_U()),
            p(_p()),
            p_rgh(_p_rgh()),
            T(_T()),
            nut(_nut()),
	    k(_k()),
	    args(_args())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& p_rgh;
        volScalarField& T;
	volScalarField& nut;
	volScalarField& k;

	/// Arg List
        argList& args;

	//void offlineSolve()
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
		ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
	        ITHACAstream::read_fields(Prghfield, p_rgh, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
		ITHACAstream::read_fields(nutFields, nut, "./ITHACAoutput/Offline/");
		//ITHACAstream::read_fields(kFields, k, "./ITHACAoutput/Offline/");
            }
            else
            {
	        for (label i = 0; i < par_BC.rows(); i++)
                {     

		     for (label j = 0; j < inletIndexGradT.size(); j++)
                     { 
		         assignBC(T, inletIndexGradT(0,0), par_BC(i, 0)); 
		     }
                      mu_now[0] = mu(0, 0);

		      truthSolve(mu_now);
                }

            }
        }


        void onlineSolveRead(fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Tfield_on, T, folder);
            }
            else
            {
            }
        }
       
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial20Steady example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./par2");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    word inlet_BC("./inlet");
    Eigen::MatrixXd inlet_Vec_BC = ITHACAstream::readMatrix(inlet_BC);
    word par_online_BC("./par3");
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesPrghproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesPrghproj", 5);
    int NmodesTproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    int NmodesNUTproj = para.ITHACAdict->lookupOrDefault<int>("NmodesNUTproj", 5);
    int NmodesKproj = para.ITHACAdict->lookupOrDefault<int>("NmodesKproj", 5);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    /// Set the parameter ranges
    example.mu_range(0, 0) = 0.000000596;
    example.mu_range(0, 1) = 0.000000596; 
    double e = para.ITHACAdict->lookupOrDefault<double>("RBFradius", 1);
    example.e = e;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    /// Set the inlet Temperature boundaries where there are non homogeneous boundary conditions
    example.inletIndexT.resize(1, 1);
    example.inletIndexT(0, 0) = 0;
    example.inletIndexGradT.resize(1, 1);
    example.inletIndexGradT(0, 0) = 4;
    // Perform the offline solve
    example.offlineSolve(par_off_BC);
    // Perform a POD decomposition for velocity temperature and pressure fields
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesOut);
    //ITHACAPOD::getModes(example.Prghfield, example.Prghmodes, example.podex, 0, 0,
    //                    NmodesOut);
    //ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                //        NmodesOut);

    //if (ITHACAutilities::check_folder("./ITHACAoutput/POD/"))
   // {
	//ITHACAPOD::getModes(example.Prghfield, example.Prghmodes, example.podex, 0, 0,
                   //     NmodesOut);
   // }
   // else 
   // {
	ITHACAPOD::getModesVelBasis(example.Prghfield, example.Prghmodes,example.Ufield, example.podex, 0, 0,
                        NmodesOut);
   // }
    
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.nutFields, example.nutModes, example.podex,
                        example.supex, 0, NmodesOut);
   // ITHACAPOD::getModes(example.kFields, example.kModes, example.podex,
              //          example.supex, 0, NmodesOut);


    // Solve the supremizer problem based on the pressure modes
    //example.solvesupremizerP("modes");
    // Solve the supremizer problem based on the shifted pressure modes
  //  example.solvesupremizerPrgh();

  //  ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
  //                      example.supex, 1,  NmodesOut);
    // Get reduced matrices
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesPrghproj, NmodesTproj,
                      NmodesSUPproj, NmodesNUTproj, NmodesKproj, par_off_BC);

    // Create the reduced object
    ReducedSteadyBBTurb reduced(example);
    // Set the inlet velocity
    Eigen::MatrixXd vel = inlet_Vec_BC;
    // Set the temperature boundary conditions
    Eigen::MatrixXd temp(1, 1);
    temp(0, 0) = 423.15;

  reduced.par_on_BC = par_on_BC.col(0);

    //reduced.rbfCoeffMat.resize(NmodesNUTproj + 1, 37);

    // Perform an online solve for the new values of inlet velocities
    for (int k = 0; k < par_off_BC.rows(); k++)
    {	
	
	reduced.gradient_t = par_off_BC.row(k);
	reduced.nu = example.mu(0, 0);
        reduced.Pr = 0.0088;
	reduced.Prt = 0.90;

        if (example.bcMethod == "penalty")
    	{
            // Set initial quess for penalty factors
	    reduced.tauIter = Eigen::MatrixXd::Zero(1,1);
            reduced.tauIter <<  1; 
	    // Solve for the penalty factors with the iterative solver
	    // Penalty factor Temperature
            reduced.tauT = reduced.tauIter*1;//1
	    // Penalty factor Velocity
	    reduced.tauU = reduced.tauIter*0.1;//0.1
	    // Penalty factor Heat Flux
	    reduced.tauHF = reduced.tauIter*5;//10
    	}
        reduced.solveOnlineSUP(vel, temp, k);
        Eigen::MatrixXd tmp_sol(reduced.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(reduced.y.rows()) = reduced.y;
        reduced.online_solution.append(tmp_sol);

	Eigen::MatrixXd tmp_sol2(NmodesNUTproj + 1, 1);
	tmp_sol2(0) = e;
	tmp_sol2.col(0).tail(reduced.rbfCoeff.rows()) = reduced.rbfCoeff;
	reduced.rbfCoeffMat.append(tmp_sol2);
    }

    // Save the online solution
    ITHACAstream::exportMatrix(reduced.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");

    ITHACAstream::exportMatrix(reduced.rbfCoeffMat, "rbfCoeffMat", "eigen",
                               "./ITHACAoutput/rbfCoeffMat");

    reduced.reconstructSUP("./ITHACAoutput/ReconstructionSUP", 1);


    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::error_listfields(
                                         example.Tfield_on, reduced.TREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");


  /*  // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::error_listfields(
                                         example.Tfield, reduced.TREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");*/
      
    return 0;
  
}

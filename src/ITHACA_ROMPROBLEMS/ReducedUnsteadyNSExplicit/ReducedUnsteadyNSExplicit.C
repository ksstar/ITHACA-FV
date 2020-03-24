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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedUnsteadyNSExplicit class


#include "ReducedUnsteadyNSExplicit.H"
#include <math.h> 

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
reducedUnsteadyNSExplicit::reducedUnsteadyNSExplicit()
{
}

reducedUnsteadyNSExplicit::reducedUnsteadyNSExplicit(unsteadyNSExplicit& FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();

    // Create locally the velocity modes
    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append(problem->supmodes[k]);
    }

    // Create locally the pressure modes
    for (label k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

}



// * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * //


void reducedUnsteadyNSExplicit::solveOnline(Eigen::MatrixXd vel,
                                        label startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");
    int numberOfStores = round(storeEvery / dt);
    vel_now = vel;

    // Create and resize the solution vectors
    Eigen::VectorXd a_n(Nphi_u);
    Eigen::VectorXd a_o(Nphi_u);

    // Initial condition
    a_o.setZero();
    a_o.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[startSnap],
                     Umodes);

    // Counting variable
    int counter = 0;

    // Set the initial time
    time = tstart;

    // Determine number of time steps
    while (time < finalTime - 0.5 * dt)
    {
        time = time + dt;
        counter ++;
    }
    time = tstart;

    // Set size of online solution
    online_solution.resize(counter+1);
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(a_o.rows()) = a_o;
    online_solution[0] = tmp_sol;

    
    for (label i = 1; i < online_solution.size(); i++)
    {
        time = time + dt;
  	std::cout << " ################## time =   " << time <<
                  " ##################" << std::endl;
    	
    	// Convective term
    	Eigen::MatrixXd cc(1, 1);
    	// Mom Term
    	Eigen::VectorXd M1 = problem->M_matrix * a_o;
	Eigen::VectorXd M2 = problem->B_matrix * a_o * nu ;

    	for (label i = 0; i < Nphi_u; i++)
    	{
        	cc = a_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_o;
         	a_n(i) =  (M2(i) - cc(0, 0) + problem->BC_matrix(i,0))*dt + a_o(i);
    	}

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(a_n.rows()) = a_n;
	online_solution[i] = tmp_sol;

	// Set oldTime for next time step.
	a_o = a_n;

    }

	ITHACAstream::exportMatrix(online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");

	int exportEveryIndex = round(exportEvery / dt);
	int outputIndex = round(finalTime / exportEvery)+1;	

}




void reducedUnsteadyNSExplicit::reconstruct(fileName folder)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }

           ITHACAstream::exportSolution(U_rec,  name(counter2), folder);
          /*   volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(P_rec, name(counter2), folder);*/
            nextwrite += exportEveryIndex;
            double timenow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append(U_rec);
           // PREC.append(P_rec);
        }

        counter++;
    }
}


//************************************************************************* //

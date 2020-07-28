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
#include "fvCFD.H"

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
    Nphi_phi = Nphi_u;

    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

    if (problem->ExplicitMethod == "A")
    {
    	for (label k = 0; k < Nphi_phi; k++)
    	{
            Phimodes.append(problem->Phimodes[k]);
    	}
    }

    newton_object = newton_unsteadyNSExplicit(Nphi_p, Nphi_p,
                        FOMproblem);

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


if (problem->ExplicitMethod == "Ales")
{
    // Create and resize the solution vectors
    Eigen::VectorXd a_n = Eigen::VectorXd::Zero(Nphi_u);
    Eigen::VectorXd a_o = Eigen::VectorXd::Zero(Nphi_u);
    Eigen::MatrixXd b_n = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::MatrixXd xx = Eigen::VectorXd::Zero(Nphi_p);

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

    a_o = ITHACAutilities::get_coeffs(problem->Ufield[0],
                     Umodes);
    b_n = ITHACAutilities::get_coeffs(problem->Pfield[0],
                     Pmodes);


    if (problem->bcMethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            a_n(j) = vel_now(j, 0);
	    a_o(j) = vel_now(j, 0);
        }
    }

    // Set size of online solution
    online_solution.resize(counter+1);
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).segment(1, Nphi_u) = a_o;
    tmp_sol.col(0).tail(b_n.rows()) = b_n;
    online_solution[0] = tmp_sol;

    /// Pressure field
    volScalarField p = problem->_p;

    for (label i = 1; i < online_solution.size(); i++)
    {

    	time = time + dt;
    	std::cout << " ################## time =   " << time <<
                  " ##################" << std::endl;

	Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(Nphi_p);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(Nphi_p);
        scalar P_norm_res(1);

	Eigen::VectorXd RHS  = Eigen::VectorXd::Zero(Nphi_p);
	// Diffusion Term
    	Eigen::VectorXd M5 = problem->DF_matrix * a_o *nu ;
	// Convection Term
	Eigen::MatrixXd cp(1, 1);
	// Mom term
	Eigen::MatrixXd M4 = problem->P_matrix * a_o;

        for (label l = 0; l < Nphi_p; l++)
    	{
      	    cp = a_o.transpose() * Eigen::SliceFromTensor(problem->Cf_tensor, 0,
                l) * a_o;
       	    RHS(l) = (1/dt) * M4(l,0)- cp(0,0) + M5(l,0);
    	}

	if (problem->bcMethod == "lift")
    	{
	    List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
	    RedLinSysP[1] =0*RedLinSysP[1];

 	    b_n = reducedProblem::solveLinearSysAxb(RedLinSysP, RHS, xx, presidual);
	}
	else
	{
	    List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
	    RedLinSysP[1] =(1/dt)*RedLinSysP[1]+nu * problem->RedLinSysPDiff[1]+ problem->RedLinSysPConv[1];

 	    b_n = reducedProblem::solveLinearSysAxb(RedLinSysP, RHS, xx, presidual);
	}


	//b_n = ITHACAutilities::get_coeffs(problem->Pfield[i],
         //            Pmodes);
	    // Convective term
    	    Eigen::MatrixXd cc(1, 1);
    	    // Diff Term
    	    Eigen::VectorXd M2 = problem->B_matrix * a_o * nu ;
	    // Pressure Term
	    Eigen::VectorXd M3 = problem->K_matrix * b_n;

	if (problem->bcMethod == "lift")
    	{

	    for (label l = 0; l < Nphi_u; l++)
    	    {
        	cc = a_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
              		l) * a_o;
        	 a_n(l) =  (M2(l)  - cc(0,0))*dt + a_o(l)- dt*M3(l);
    	    }
            for (label j = 0; j < N_BC; j++)
            {
            	a_n(j) = vel_now(j, 0);
            }
    	}
	else
	{
            for (label l = 0; l < Nphi_u; l++)
    	    {
        	cc = a_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
              		l) * a_o;
        	 a_n(l) =  (M2(l)  - cc(0,0))*dt  +problem->BC_matrix(l,0)*dt+ a_o(l)- dt*M3(l);
    	    }
	}


	//a_n = ITHACAutilities::get_coeffs(problem->Ufield[i],
        //             Umodes); 
    	tmp_sol(0) = time;
    	tmp_sol.col(0).segment(1, Nphi_u) = a_n;
    	tmp_sol.col(0).tail(b_n.rows()) = b_n;
    	online_solution[i] = tmp_sol;
 
    	a_o = a_n;
    }
}
else if (problem->ExplicitMethod == "A")
{

    // Create and resize the solution vectors
    Eigen::VectorXd a_n = Eigen::VectorXd::Zero(Nphi_u);
    Eigen::VectorXd a_o = Eigen::VectorXd::Zero(Nphi_u);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::MatrixXd xx = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::VectorXd c_o = Eigen::VectorXd::Zero(Nphi_phi);
    Eigen::VectorXd c_n = Eigen::VectorXd::Zero(Nphi_phi);


    a_o = ITHACAutilities::get_coeffs(problem->Ufield[0],
                   Umodes);
    b = ITHACAutilities::get_coeffs(problem->Pfield[0],
                   Pmodes);
    c_o = ITHACAutilities::get_coeffs(problem->Phifield[0],
                   Phimodes);


     if (problem->bcMethod == "lift")
    	{
            for (label j = 0; j < N_BC; j++)
            {
            	a_n(j) = vel_now(j, 0);
	   	a_o(j) = vel_now(j, 0);
            }
    	}

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
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + Nphi_u + 1, 1);

    tmp_sol(0) = time;
    tmp_sol.col(0).segment(1, Nphi_u) = a_o;
    tmp_sol.col(0).segment(Nphi_u + 1, Nphi_p) = b;
    tmp_sol.col(0).tail(Nphi_u) = c_o;
    online_solution[0] = tmp_sol;


 //   tmp_sol(0) = time;
  //  tmp_sol.col(0).segment(1, Nphi_u) = a_o;
  //  tmp_sol.col(0).tail(b.rows()) = b;
 //   online_solution[0] = tmp_sol;


    for (label i = 1; i < online_solution.size(); i++)
    {
	Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(Nphi_p);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(Nphi_p);
        scalar P_norm_res(1);
        time = time + dt;
  	std::cout << " ################## time =   " << time <<
                  " ##################" << std::endl;

	Eigen::VectorXd RHS  = Eigen::VectorXd::Zero(Nphi_p);
	// Diffusion Term
    	Eigen::VectorXd M5 = problem->DF_matrix * a_o *nu ;
	// Convection Term
	Eigen::MatrixXd cp(1, 1);
	// Mom term
	Eigen::MatrixXd M4 = problem->P_matrix * a_o;

    	for (label l = 0; l < Nphi_p; l++)
    	{
        	cp = c_o.transpose() * Eigen::SliceFromTensor(problem->Cf_tensor, 0,
                l) * a_o;
        	RHS(l) = (1/dt) * M4(l,0) - cp(0,0)+ M5(l,0);

    	}
	//List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
	//RedLinSysP[1] =(1/dt)*RedLinSysP[1]+ problem->RedLinSysPConv[1]+ nu * problem->RedLinSysPDiff[1];
	//b = reducedProblem::solveLinearSysAxb(RedLinSysP, RHS, xx, presidual);
    	//b = reducedProblem::solveLinearSysAxb(problem->RedLinSysP, RHS, xx, presidual);

	if (problem->bcMethod == "lift")
    	{
		List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
	   	RedLinSysP[1] =0*RedLinSysP[1];

 	   	b = reducedProblem::solveLinearSysAxb(RedLinSysP, RHS, xx, presidual);
	}
    	else
	{
		List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
	   	RedLinSysP[1] =(1/dt)*RedLinSysP[1]+nu * problem->RedLinSysPDiff[1]+ problem->RedLinSysPConv[1];

 	   	b = reducedProblem::solveLinearSysAxb(RedLinSysP, RHS, xx, presidual);
	}

	// Convective term
    	Eigen::MatrixXd cc(1, 1);
    	// Diff Term
    	Eigen::VectorXd M2 = problem->B_matrix * a_o * nu ;
	// Pressure Term
	Eigen::VectorXd M3 = problem->K_matrix * b;

	 if (problem->bcMethod == "lift")
    	{
	    for (label l = 0; l < Nphi_u; l++)
    	    {
        	cc = c_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
              		l) * a_o;
        	 a_n(l) = a_o(l) + (M2(l)  - cc(0,0)- M3(l))*dt  ;
    	    }
            for (label j = 0; j < N_BC; j++)
            {
            	a_n(j) = vel_now(j, 0);
            }
    	}
	else
	{
	    for (label l = 0; l < Nphi_u; l++)
    	    {
        	cc = c_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
              		l) * a_o;
        	 a_n(l) = a_o(l) + (M2(l)  - cc(0,0)+problem->BC_matrix(l,0)- M3(l))*dt  ;
    	    }
	}

	Eigen::MatrixXd M6 = problem->I_matrix * a_o;
	Eigen::MatrixXd M7 = problem->ID_matrix * a_o*nu;
	Eigen::MatrixXd M8 = problem->Kf_matrix*b.col(0);
	Eigen::MatrixXd M11(Nphi_phi,Nphi_phi);
	Eigen::MatrixXd M12 = Eigen::VectorXd::Zero(Nphi_phi);
	Eigen::MatrixXd ci(Nphi_phi,1);
	
	for (label k = 0; k < Nphi_phi; k++)
    	{	
		M12 = dt*problem->Ci_matrix[k] * a_o*c_o(k) +M12;//
	}

	if (problem->bcMethod == "lift")
    	{
		c_n = problem->Mf_matrix.colPivHouseholderQr().solve(M6-M12+ dt*(-M8+M7));
	}
	else
	{
		c_n = problem->Mf_matrix.colPivHouseholderQr().solve(M6-M12+ dt*(-M8+M7
					+nu*problem->BC_matrix_PPE-problem->Bbc_matrix));
	}

	tmp_sol(0) = time;
    	tmp_sol.col(0).segment(1, Nphi_u) = a_n;
	tmp_sol.col(0).segment(Nphi_u + 1, Nphi_p) = b;
    	tmp_sol.col(0).tail(Nphi_u) = c_n;
    	online_solution[i] = tmp_sol;

	a_o = a_n;
	c_o = c_n;


    	}
}

	/*ITHACAstream::exportMatrix(online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");

	int exportEveryIndex = round(exportEvery / dt);
	int outputIndex = round(finalTime / exportEvery)+1;	*/
	
}




void reducedUnsteadyNSExplicit::reconstruct(fileName folder)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);

    Eigen::MatrixXd ContErrROM;
    ContErrROM.resize(online_solution.size(),2);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
	
	  //  ITHACAstream::exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

	  // ITHACAstream::exportSolution(P_rec, name(counter2), folder);

	    if (problem->ExplicitMethod == "Ales")
    	    {

	       volScalarField contErr = fvc::div(U_rec);
   	       ContErrROM(i,0) = dt*
        	   mag(contErr)().weightedAverage(contErr.mesh().V()).value(); //scalar sumLocalContErr
	       ContErrROM(i,1) = dt*
        	   contErr.weightedAverage(contErr.mesh().V()).value(); //scalar globalContErr
	   }
	   else if (problem->ExplicitMethod == "A")
    	    {
	    	surfaceScalarField Phi_rec("Phi_rec", Phimodes[0] * 0);

            	for (label j = 0; j < Nphi_u; j++)
            	{
                    Phi_rec += Phimodes[j] * online_solution[i](j + Nphi_u + Nphi_p + 1, 0);
            	}

	   volScalarField contErr = fvc::div(Phi_rec);
   	   ContErrROM(i,0) = dt*
        	mag(contErr)().weightedAverage(contErr.mesh().V()).value(); //scalar sumLocalContErr
	   ContErrROM(i,1) = dt*
        	contErr.weightedAverage(contErr.mesh().V()).value(); //scalar globalContErr

           // ITHACAstream::exportSolution(Phi_rec,  name(counter2), folder);
	   }


            nextwrite += exportEveryIndex;
            double timenow = online_solution[i](0, 0);
           // std::ofstream of(folder + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }

	ITHACAstream::exportMatrix(ContErrROM, "ContErrROM" + name(Nphi_u), "eigen", folder);
}



//************************************************************************* //

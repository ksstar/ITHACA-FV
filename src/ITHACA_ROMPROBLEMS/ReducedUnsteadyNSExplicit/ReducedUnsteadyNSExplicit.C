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

    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

    newton_object = newton_unsteadyNSExplicit(Nphi_p, Nphi_p,
                        FOMproblem);

}

// * * * * * * * * * * * *  Operators Poisson FOM level   * * * * * * * * * * //

// Operator to evaluate the residual for the Pressure Poisson Equation (PPE) approach
	int newton_unsteadyNSExplicit::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
	{

	Eigen::VectorXd b_tmp(Nphi_p);
        b_tmp = x;
	//Eigen::MatrixXd M5 = (1/dt)*problem->PF_matrix * a_n;
	//Eigen::MatrixXd M6 = (1/dt)*problem->DF_matrix * a_n;

 	// Pressure Term
    	Eigen::VectorXd M5 = problem->D_matrix * b_tmp;

    	// BC PPE
    	Eigen::VectorXd M6 = problem->BC3_matrix * a_n * nu;

	Eigen::VectorXd M7 = (problem->K_matrix).transpose() * a_n;

	Eigen::MatrixXd gg(1, 1);

    	for (label l = 0; l < Nphi_p; l++)
    	{
        	gg = a_n.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0,
                l) * a_n;
        	fvec(l) = -(1/dt) * M7(l,0) + M5(l, 0) + gg(0, 0) - M6(l, 0);

    	}

    	return 0;
	}

	// Operator to evaluate the Jacobian for the supremizer approach
	int newton_unsteadyNSExplicit::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
	{
    	Eigen::NumericalDiff<newton_unsteadyNSExplicit> numDiff(*this);
    	numDiff.df(x, fjac);
    	return 0;
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
  

    Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::VectorXd presidual = Eigen::VectorXd::Zero(Nphi_p);
    scalar P_norm_res(1);


   /* y.resize(Nphi_p, 1);
    y.setZero();
    newton_object.dt = dt;
    newton_object.nu = nu;
    newton_object.a_n = a_n;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNSExplicit> hnls(newton_object);*/
    
    

    for (label i = 1; i < online_solution.size(); i++)
    {

    	time = time + dt;
    	std::cout << " ################## time =   " << time <<
                  " ##################" << std::endl;


	if (problem->PoissonMethod == "ROM")
	{
    	 /*  Eigen::MatrixXd M4 = (1/dt)*(-1)*(problem->K_matrix).transpose() * a_n;
    	   b_n = reducedProblem::solveLinearSysAxb(problem->RedLinSysP, M4, xx, presidual);

    	   Eigen::VectorXd M3 = problem->K_matrix * b_n;

    	   for (label k = 0; k < Nphi_u; k++)
    	   {  
		a_n(k) =  a_n(k) - (dt * M3(k));
    	   }*/
	}
	else if (problem->PoissonMethod == "FOM")
	{

	Eigen::VectorXd RHS  = Eigen::VectorXd::Zero(Nphi_p);
	// Diffusion Term
    	Eigen::VectorXd M5 = problem->DF_matrix * a_o * nu;
	// Convection Term
	Eigen::MatrixXd cp(1, 1);
	// Mom term
	Eigen::MatrixXd M4 = problem->P_matrix * a_o;

    	for (label l = 0; l < Nphi_p; l++)
    	{
        	cp = a_o.transpose() * Eigen::SliceFromTensor(problem->Cf_tensor, 0,
                l) * a_o;
        	RHS(l) = (1/dt) * M4(l,0)  -cp(0,0)+ M5(l,0);

    	}

    	b_n = reducedProblem::solveLinearSysAxb(problem->RedLinSysP, RHS, xx, presidual);

	// Convective term
    	Eigen::MatrixXd cc(1, 1);
    	// Mom Term
    	Eigen::VectorXd M1 = problem->M_matrix * a_o;
    	// Diff Term
    	Eigen::VectorXd M2 = problem->B_matrix * a_o * nu ;
	// Pressure Term
	Eigen::VectorXd M3 = problem->K_matrix * b_n;

    	for (label l = 0; l < Nphi_u; l++)
    	{
        	cc = a_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
              		l) * a_o;
        	 a_n(l) =  (M2(l)  - cc(0,0))*dt -problem->BC_matrix(l,0)*dt + a_o(l) - dt*M3(l);
    	}

    	  /* 

    	   for (label k = 0; k < Nphi_u; k++)
    	   {  
		a_n(k) =  a_n(k) - (dt * M3(k));
    	   }*/
	 

	   /*newton_object.nu = nu;
	   newton_object.a_n = a_n;

	   newton_object.dt = dt;
	   
	    Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);

	    newton_object.operator()(y, res);
	    b_n = y;*/

	}

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
    Eigen::MatrixXd x = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::VectorXd c_o = Eigen::VectorXd::Zero(Nphi_phi);
    Eigen::VectorXd c_n = Eigen::VectorXd::Zero(Nphi_phi);


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
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).segment(1, Nphi_u) = a_o;
    tmp_sol.col(0).tail(b.rows()) = b;
    online_solution[0] = tmp_sol;

    /// Pressure field
    volScalarField& p = problem->_p();

    Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(Nphi_p);
    Eigen::VectorXd presidual = Eigen::VectorXd::Zero(Nphi_p);
    scalar P_norm_res(1);

    for (label i = 1; i < online_solution.size(); i++)
    {
	presidual = presidual*0;
        time = time + dt;
  	std::cout << " ################## time =   " << time <<
                  " ##################" << std::endl;

	
    	// Convective term
    	Eigen::MatrixXd cc(1, 1);
    	// Mom Term
    	Eigen::VectorXd M1 = problem->M_matrix * a_o;
	Eigen::VectorXd M2 = problem->B_matrix * a_o * nu ;

    	for (label l = 0; l < Nphi_u; l++)
    	{
        	cc = c_o.transpose() * Eigen::SliceFromTensor(problem->Cf_tensor, 0,
                l) * a_o;
         	a_n(l) =  (M2(l)  - cc(0,0)-problem->BC_matrix(l,0))*dt + a_o(l);
    	}

	
	Eigen::VectorXd M6 = problem->I_matrix * a_n;

	for (label k = 0; k < Nphi_p; k++)
    	{
		c_n(k) =  M6(k) ;
	}

	Eigen::MatrixXd M4 = (1/dt)*problem->Pf_matrix * a_n;
        b = reducedProblem::solveLinearSysAxb(problem->RedLinSysP, M4, x, presidual);

	Eigen::VectorXd M3 = problem->K_matrix * b;

	for (label k = 0; k < Nphi_u; k++)
    	{
		a_n(k) =  a_n(k) - (dt * M3(k));
	}

	Eigen::VectorXd M5 = problem->Kf_matrix * b;
	
	for (label k = 0; k < Nphi_phi; k++)
    	{
		c_n(k) =  M6(k) - (dt * M5(k));
	}

	tmp_sol(0) = time;
    	tmp_sol.col(0).segment(1, Nphi_u) = a_n;
    	tmp_sol.col(0).tail(b.rows()) = b;
    	online_solution[i] = tmp_sol;
 
	a_o = a_n;
	c_o = c_n;

    	}
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

	forAll( U_rec.boundaryFieldRef()[0], l)
	{
		U_rec.boundaryFieldRef()[0][l].component(vector::X) = 1;
	}

           ITHACAstream::exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

	   // volVectorField P_grad= fvc::grad(P_rec);
            ITHACAstream::exportSolution(P_rec, name(counter2), folder);
	   //ITHACAstream::exportSolution(P_grad, name(counter2), folder);
            nextwrite += exportEveryIndex;
            double timenow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }
}


//************************************************************************* //

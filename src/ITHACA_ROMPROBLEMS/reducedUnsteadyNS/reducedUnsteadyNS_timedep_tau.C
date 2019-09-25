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
/// Source file of the reducedUnsteadyNS class


#include "reducedUnsteadyNS.H"

#include <cmath> 
#include <iomanip> 
#include <iostream> 
  
using namespace std; 


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
reducedUnsteadyNS::reducedUnsteadyNS()
{
}

reducedUnsteadyNS::reducedUnsteadyNS(unsteadyNS& FOMproblem)
{
    problem = &FOMproblem;
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();

    // Create locally the velocity modes
    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        LUmodes.append(problem->liftfield[k]);
    }

    if (problem->aveMethod == "mean")
    {
       //scalar magSumSquare = Foam::sqrt(fvc::domainIntegrate(problem->Umean[0] & problem->Umean[0]).value());
       //LUmodes.append((problem->Umean[0]/magSumSquare));
       LUmodes.append(problem->Umean[0]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        LUmodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        LUmodes.append(problem->supmodes[k]);
    }

    if (problem->aveMethod == "mean")
    {
       MPmodes.append(problem->Pmean[0]);
    }

    for (label k = 0; k < problem->NPmodes; k++)
    {
       MPmodes.append(problem->Pmodes[k]);
    }

    newton_object_sup = newton_unsteadyNS_sup(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
    newton_object_PPE = newton_unsteadyNS_PPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
}

// * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Supremizer approach
int newton_unsteadyNS_sup::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (label l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l,0) * BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              a_tmp;
        }
    }

    for (label i = 0; i < Nphi_u; i++)
    {
   	cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);

	if (problem->bcMethod == "penalty")
        {
            fvec(i) += (penaltyU.row(i).col(0))(0, 0);
	}
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_sup::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Pressure Poisson Equation (PPE) approach
int newton_unsteadyNS_PPE::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{

    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;

    // Split K_matrix using Gauss
   // Eigen::VectorXd M2A = problem->K1_matrix * b_tmp;
   // Eigen::VectorXd M2B = problem->K2_matrix * b_tmp;

    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;
    // BC PPE
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;
    // BC PPE - penalty method
    Eigen::VectorXd M8 = problem->PF_matrix * a_dot;


    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (label l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = abs(tauU(l,0)) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              a_tmp);
        }
    }

    for (label i = 0; i < Nphi_u; i++)
    {

      cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;

        //fvec(i) = - M5(i) + M1(i) - cc(0, 0) + M2A(i) - M2B(i) ;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i) ;

  	if (problem->bcMethod == "penalty")
        {
	    for (label l = 0; l < N_BC; l++)
            {
                 fvec(i) += penaltyU(i,l);
	    }
	}
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        gg = a_tmp.transpose() * problem->G_matrix[j] * a_tmp;

        fvec(k) = M3(j, 0) + gg(0, 0) - M7(j, 0)+ M8(j,0);
	
    }

    if (problem->bcMethod == "lift")
    {
    	for (label j = 0; j < N_BC; j++)
    	{
            fvec(j) = x(j) - BC(j);
    	}
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_PPE::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * Solve Functions supremizer * * * * * * * * * * * //

Eigen::MatrixXd reducedUnsteadyNS::solveOnline_sup(Eigen::MatrixXd& vel_now, label NParaSet,
                                        label startSnap)
{

    std::cout << "################## Online solve N° " << NParaSet <<
    " ##################" << std::endl;
    std::cout << "Solving for the parameter: " << vel_now << std::endl;
    // Count number of time steps
    int counter = 0;
    time = tstart;
    while (time < finalTime - 0.5 * dt)
    {
        time = time + dt;
        counter ++;
    }
    // Set size of online solution
    online_sol.resize(Nphi_u + Nphi_p + 1, counter + 1);
    // Set initial condition for online solve
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p,  1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0], LUmodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[0], MPmodes);

    // Set some properties of the newton objects
    newton_object_sup.nu = nu;
    newton_object_sup.y_old = y;
    newton_object_sup.dt = dt;
    newton_object_sup.BC.resize(N_BC);
    newton_object_sup.tauU = tauU;

    for (label j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }

    // Set the initial time
    time = tstart;

    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_sol.col(0) = tmp_sol;

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

  // Start the time loop
    for (label i = 1; i < online_sol.cols(); i++)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        if (problem->bcMethod == "lift")
        {
            for (label j = 0; j < N_BC; j++)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old = y;
        Info << "Time = " << time << endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
            hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
            hnls.iter << " iterations " << def << std::endl << std::endl;
        }

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        online_sol.col(i) = tmp_sol;
    }

    // Save the current solution
    ITHACAstream::exportMatrix(online_sol, "red_coeff", "python",
       "./ITHACAoutput/red_coeff/" + name(NParaSet) + "/");
   // ITHACAstream::exportMatrix(online_sol, "red_coeff", "matlab",
     //  "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
   // ITHACAstream::exportMatrix(online_sol, "red_coeff", "eigen",
     //  "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    return online_sol;
}


// * * * * * * * * * * * * * * * Solve Functions PPE * * * * * * * * * * * * * //

Eigen::MatrixXd reducedUnsteadyNS::solveOnline_PPE(Eigen::MatrixXd& vel, label NParaSet,
                                        label startSnap)
{
    std::cout << "################## Online solve N° " << NParaSet <<
    " ##################" << std::endl;
    std::cout << "Solving for the parameter: " << vel.col(0) << std::endl;

    if (problem->bcMethod == "lift")
    {
 	vel_now = setOnlineVelocity(vel);
    }
    else
    {
   	vel_now = vel;
    }

     // Count number of time steps
    int counter = 0;
    time = tstart;

    while (time < finalTime - 0.5 * dt)
    {
        time = time + dt;
        counter ++;
    }

    // Set size of online solution
    online_sol.resize(Nphi_u + Nphi_p + 1, counter + 1);

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0], LUmodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[0], MPmodes);
   
    if (problem->bcMethod == "lift")
    {
    	// Change initial condition for the lifting function
    	for (label j = 0; j < N_BC; j++)
    	{
        	y(j) = vel_now(j, 0);
     	}
    }

    tauIter = tauInit;

    // Set some properties of the newton object
    newton_object_PPE.nu = nu;
    newton_object_PPE.y_old = y;
    newton_object_PPE.dt = dt;
    newton_object_PPE.BC.resize(N_BC);
    newton_object_PPE.tauU = tauIter;

    for (label j = 0; j < N_BC; j++)
    {
            newton_object_PPE.BC(j) = vel_now(j, 0);
    }

    // Set the initial time
    time = tstart;
    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    online_sol.col(0) = tmp_sol;
  
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(N_BC,counter + 1);
    // initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(N_BC,counter + 1);

      // Start the time loop
    for (label i = 1; i < online_sol.cols(); i++)
    {

      time = time + dt;

  //  std::cout << "Beginning time loop" << std::endl;
    int Iter = 0;

    Eigen::VectorXd diffvel = vel_now.col(i) - valBC.col(i);

 //   std::cout << "diffvel(0,0)" << diffvel(0,0) << std::endl;

    	while (diffvel.maxCoeff()> tolerance && Iter < maxIter)
    	{

 //std::cout << "Beginning while loop" << std::endl;
            if (Iter == 0)
            {
                tauIter = tauIter; 
            }
           else
           {

                for (label j = 0; j < N_BC; j++)
                {
		     tauIter(j,0) = tauIter(j,0)*diffvel(j)/(tolerance);

	        }

           }

	//std::cout << "tauIter(0,0): " << tauIter(0,0) << std::endl;

	   newton_object_PPE.tauU = tauIter;

    	   //  Set the old boundary value to the current value
           valBC0.col(i)  = valBC.col(i);

           if (problem->timedepbcMethod == "yes" )
           {
	    for (label j = 0; j < N_BC; j++)
            {
                newton_object_PPE.BC(j) = vel_now(j, i);
            }
           }

           Eigen::VectorXd res(y);
           res.setZero();
           hnls.solve(y);
           if (problem->bcMethod == "lift")
           {

	    if (problem->timedepbcMethod == "yes" )
            {
  		for (label j = 0; j < N_BC; j++)
                {
                    y(j) = vel_now(j, i);
                }
	    }

	    else if (problem->timedepbcMethod == "no" )
	    {
                for (label j = 0; j < N_BC; j++)
                {
                    y(j) = vel_now(j, 0);
                }
	    } 
           }
           newton_object_PPE.operator()(y, res);
           newton_object_PPE.y_old = y;
           Info << "Time = " << time << endl;

          if (res.norm() < 1e-5)
          {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
            hnls.iter << " iterations " << def << std::endl << std::endl;
          }
          else
          {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
            hnls.iter << " iterations " << def << std::endl << std::endl;
          }

	  volVectorField U_rec("U_rec", LUmodes[0] * 0);
          for (label j = 0; j < Nphi_u; j++)
          {
                U_rec += LUmodes[j] * y(j);
          }

          for (label k = 0; k < problem->inletIndex.rows(); k++)
          {
                 label BCind = problem->inletIndex(k,0);
		 label BCcomp = problem->inletIndex(k,1);
                 valBC(k,i) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);
		 diffvel(k) = abs(abs(vel_now(k,i)) - abs(valBC(k,i)));
	  }

	   //std::cout << "max error: " << diffvel.maxCoeff() << std::endl;

          // Count the number of iterations
          Iter ++;

       } // while loop

    std::cout << "Final penalty factor(s): " << tauIter << std::endl;
    std::cout << "Iterations: " << Iter << std::endl;

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        online_sol.col(i) = tmp_sol;
    }

    ITHACAstream::exportMatrix(online_sol, "red_coeff", "eigen",
       "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    return online_sol;
 
}


// * * * * * * * * * * * * * * * Calculate penalty factor * * * * * * * * * * * * * //
Eigen::MatrixXd reducedUnsteadyNS::penalty_sup(Eigen::MatrixXd& vel_now, Eigen::MatrixXd& tauInit)
{
    // initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(vel_now.rows(),1);
    // initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(vel_now.rows(),1);

    tauIter = tauInit;         
    int Iter = 0;

    while (abs((vel_now - valBC).sum()) > tolerance && Iter < maxIter)
    {
        std::cout << "diff: " << abs((vel_now - valBC).sum()) << std::endl;

        if ((valBC - valBC0).sum() != 0)
        {
            for (label j = 0; j < N_BC; j++)
            {

                Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(Nphi_u, 1);

                Jacobian = vel_now(j,0)* problem->bcVelVec[j] - 
                    problem->bcVelMat[j]* y.head(Nphi_u) ;

                std::cout << "Jacobian: " << Jacobian(0,0) << std::endl;

                tauIter(j,0) = tauIter(j,0) - (valBC(j,0) - vel_now(j,0))/Jacobian(0,0);
            }
        }
        else
        {  
        }

        std::cout << "Solving for penalty factor(s): " << tauIter << std::endl;

    //  Set the old boundary value to the current value
        valBC0  = valBC;
    // Count the number of iterations
        Iter ++;

        y.resize(Nphi_u + Nphi_p, 1);
        y.setZero();
        y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0],LUmodes);
        y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[0],MPmodes);

    // Set some properties of the newton object
        newton_object_sup.nu = nu;
        newton_object_sup.y_old = y;
        newton_object_sup.dt = dt;
        newton_object_sup.BC.resize(N_BC);
        newton_object_sup.tauU = tauIter;

    // Change initial condition for the lifting function
        newton_object_sup.BC(0) = vel_now(0, 0);
 
    // Create nonlinear solver object
        Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
        Color::Modifier red(Color::FG_RED);
        Color::Modifier green(Color::FG_GREEN);
        Color::Modifier def(Color::FG_DEFAULT);

    // Set initially for convergence check
        Eigen::VectorXd res(y);
        res.setZero();

    // Start the time loop
        for (label i = 1; i < timeSteps; i++)
        {
            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);

            newton_object_sup.operator()(y, res);
            newton_object_sup.y_old = y;

            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            volVectorField U_rec("U_rec", LUmodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * y(j);
            }


            for (label k = 0; k < problem->inletIndex.rows(); k++)
            {
                 label BCind = problem->inletIndex(k,0);
		 label BCcomp = problem->inletIndex(k,1);
                 valBC(0,0) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);
            }
            std::cout << "valBC: "<< valBC << std::endl; 
        } 

        // Check whether solution has been converged
        if (res.norm() > 1e-5)
        {
            std::cout << "Penalty method has failed" << std::endl;
            exit(0);
        }
        
    }

    std::cout  << "tauIter: "<< tauIter << std::endl;


    return tauIter;
}




Eigen::MatrixXd reducedUnsteadyNS::penalty_PPE(Eigen::MatrixXd& vel_now, Eigen::MatrixXd& tauInit)
{
    // initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(N_BC,1);
    // initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(N_BC,1);

    tauIter = tauInit;         
    int Iter = 1;
    Eigen::VectorXd diffvel = vel_now.col(0) - valBC.col(0);

    while (diffvel.maxCoeff()> tolerance && Iter < maxIter)
    {

        if ((valBC - valBC0).sum() == 0)
        {
            tauIter = tauIter; 
        }
        else
        {

            for (label j = 0; j < N_BC; j++)
            {
		//tauIter(j,0) = -0.5*tauIter(j,0) + 0.5*std::sqrt(abs(tauIter(j,0)*tauIter(j,0) +
				//4*abs(valBC(j,0) - vel_now(j,0))));
		//ltauIter(j,0) = tauIter(j,0)*Iter;
		tauIter(j,0) = tauIter(j,0)*diffvel(j)/tolerance;
            }
        }

        std::cout << "Solving for penalty factor(s): " << tauIter << std::endl;

    //  Set the old boundary value to the current value
        valBC0  = valBC;

        y.resize(Nphi_u + Nphi_p, 1);
        y.setZero();
        y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0],LUmodes);
        y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[0],MPmodes);
      
    // Set some properties of the newton object
        newton_object_PPE.nu = nu;
        newton_object_PPE.y_old = y;
        newton_object_PPE.dt = dt;
        newton_object_PPE.BC.resize(N_BC);
        newton_object_PPE.tauU = tauIter;

    // Set boundary conditions
        for (label j = 0; j < N_BC; j++)
        {
            newton_object_PPE.BC(j) = vel_now(j, 0);
        }
     
    // Create nonlinear solver object
        Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
        Color::Modifier red(Color::FG_RED);
        Color::Modifier green(Color::FG_GREEN);
        Color::Modifier def(Color::FG_DEFAULT);

    // Set initially for convergence check
        Eigen::VectorXd res(y);
        res.setZero();

    // Start the time loop
        for (label i = 1; i < timeSteps; i++)
        {

	    // Set boundary conditions
            for (label j = 0; j < N_BC; j++)
            {
                newton_object_PPE.BC(j) = vel_now(j, 0);
            }

            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);

            newton_object_PPE.operator()(y, res);
            newton_object_PPE.y_old = y;

            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            volVectorField U_rec("U_rec", LUmodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * y(j);
            }

            for (label k = 0; k < problem->inletIndex.rows(); k++)
            {
                 label BCind = problem->inletIndex(k,0);
		 label BCcomp = problem->inletIndex(k,1);
                 valBC(k,0) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);
            }
            
        } 

           std::cout << "valBC: "<< valBC << std::endl;

	    for (label j = 0; j < N_BC; j++)
            {
		  diffvel(j) = abs(abs(vel_now(j,0)) - abs(valBC(j,0)));
            }


	std::cout << "max tolerance: " << diffvel << std::endl;
        std::cout << "Iterations: " << Iter << std::endl;


        // Count the number of iterations
        Iter ++;


        // Check whether solution has been converged
        if (res.norm() > 1e-5)
        {
            std::cout << "Penalty method has failed" << std::endl;
            exit(0);
        }

        
    }

std::cout << "Final penalty factor(s): " << tauIter << std::endl;
std::cout << "Iterations: " << Iter-1 << std::endl;

    return tauIter;
}




Eigen::MatrixXd reducedUnsteadyNS::penalty_PPE_time(Eigen::MatrixXd& vel_now, Eigen::MatrixXd& tauInit)
{

    // initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(N_BC,timeSteps);
    // initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(N_BC,timeSteps);

    tauIter = tauInit;         
    int Iter = 0;

    Eigen::VectorXd diffvel = vel_now.col(timeSteps-1) - valBC.col(timeSteps-1);
 std::cout << "diffvel =  " << diffvel << std::endl;
//std::cout << "diffvel.norm() =  " << diffvel.norm() << std::endl;
    while (diffvel.maxCoeff()> tolerance && Iter < maxIter)
// while (abs((vel_now.col(timeSteps-1) - valBC.col(timeSteps-1))).maxCoeff() > tolerance && Iter < maxIter)
    {
        if ((valBC.col(timeSteps-1) - valBC0.col(timeSteps-1)).sum() == 0)
        {
            tauIter = tauIter; 
        }
        else
        {

            for (label j = 0; j < N_BC; j++)
            {

		//tauIter(j,0) = 0.5*tauIter(j,0) + std::sqrt(abs(tauIter(j,0)*tauIter(j,0) +
		//		4*abs(valBC(j,0) - vel_now(j,0))));

	        //tauIter(j,0) = 2*tauIter(j,0);

		tauIter(j,0) = tauIter(j,0)*diffvel(j)/tolerance;



	    }

        }

        std::cout << "Solving for penalty factor(s): " << tauIter << std::endl;
        std::cout << "number of iterations: " << Iter << std::endl;

    //  Set the old boundary value to the current value
        valBC0  = valBC;

        y.resize(Nphi_u + Nphi_p, 1);
        y.setZero();
        y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0],LUmodes);
        y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[0],MPmodes);
 
    // Set some properties of the newton object
        newton_object_PPE.nu = nu;
        newton_object_PPE.y_old = y;
        newton_object_PPE.dt = dt;
        newton_object_PPE.BC.resize(N_BC);
        newton_object_PPE.tauU = tauIter;

    // Set boundary conditions
        for (label j = 0; j < N_BC; j++)
        {
            newton_object_PPE.BC(j) = vel_now(j, 0);
        }

    // Create nonlinear solver object
        Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
        Color::Modifier red(Color::FG_RED);
        Color::Modifier green(Color::FG_GREEN);
        Color::Modifier def(Color::FG_DEFAULT);

    // Set initially for convergence check
        Eigen::VectorXd res(y);
        res.setZero();

    // Start the time loop
        for (label i = 1; i < timeSteps; i++)
        {

	    // Set boundary conditions
            for (label j = 0; j < N_BC; j++)
            {
                newton_object_PPE.BC(j) = vel_now(j, i);
            }

            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);

            newton_object_PPE.operator()(y, res);
            newton_object_PPE.y_old = y;

            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                hnls.iter << " iterations " << def << std::endl << std::endl;
            }

            volVectorField U_rec("U_rec", LUmodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * y(j);
            }

            for (label k = 0; k < problem->inletIndex.rows(); k++)
            {
                 label BCind = problem->inletIndex(k,0);
		 label BCcomp = problem->inletIndex(k,1);
                 valBC(k,i) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);

            }
           
        } 


        // Check whether solution has been converged
        if (res.norm() > 1e-5)
        {
            std::cout << "Penalty method has failed" << std::endl;
            exit(0);
        }

    	for (label j = 0; j < N_BC; j++)
    	{
		diffvel(j) = abs(abs(vel_now(j,timeSteps-1)) - abs(valBC(j,timeSteps-1)));
    	}

	std::cout << "max error: " << diffvel.maxCoeff() << std::endl;

        // Count the number of iterations
        Iter ++;
        
    }

std::cout << "Final penalty factor(s): " << tauIter << std::endl;
std::cout << "Iterations: " << Iter-1 << std::endl;

    return tauIter;
}

// * * * * * * * * * * * * * * * Reconstruct fields * * * * * * * * * * * * * //
void reducedUnsteadyNS::reconstruct_PPE(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;

int counter2 = 1 + UREC.size();

    for (label i = 0; i < online_sol.cols(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", LUmodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * online_sol(j + 1, i);
            }

            ITHACAstream::exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", MPmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                    P_rec += MPmodes[j] * online_sol(j + Nphi_u + 1, i);
            }

            ITHACAstream::exportSolution(P_rec,  name(counter2), folder);
           
            nextwrite += printevery;

            double timenow = online_sol(0, i);

            std::ofstream of(folder + "/" + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }

}

void reducedUnsteadyNS::reconstruct_sup(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1 + UREC.size();

for (label i = 0; i < online_sol.cols(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", LUmodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * online_sol(j + 1, i);
            }

            ITHACAstream::exportSolution(U_rec,  name(counter2), folder);

	//    surfaceScalarField Phi_rec("Phi_rec", (linearInterpolate(LUmodes[0] ) & LUmodes[0].mesh().Sf() )* 0);
	//    Phi_rec = linearInterpolate(U_rec ) & U_rec.mesh().Sf();

	//    ITHACAstream::exportSolution(Phi_rec, name(counter2), folder);

 
            volScalarField P_rec("P_rec", MPmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += MPmodes[j] * online_sol(j + Nphi_u + 1, i);
            }

            ITHACAstream::exportSolution(P_rec,  name(counter2), folder);
      
       
            nextwrite += printevery;


            double timenow = online_sol(0, i);

            std::ofstream of(folder + "/" + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append(U_rec);
	    PREC.append(P_rec);
        }

        counter++;
    }

}


Eigen::MatrixXd reducedUnsteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;

    vel_scal.resize(vel.rows(), vel.cols());

    for (int i = 0; i < vel.cols(); i++)
    {

    	for (int k = 0; k < problem->inletIndex.rows(); k++)
    	{
            label p = problem->inletIndex(k, 0);
            label l = problem->inletIndex(k, 1);
            scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
            scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;
            vel_scal(k, i) = vel(k, i) / abs(u_lf);
    	}
    }

    return vel_scal;
}

// ************************************************************************* //

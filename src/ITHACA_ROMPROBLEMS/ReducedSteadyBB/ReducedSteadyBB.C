/*---------------------------------------------------------------------------*\
v/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2019 by the ITHACA-FV authors
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
/// Source file of the ReducedSteadyBB class


#include "ReducedSteadyBB.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedSteadyBB::ReducedSteadyBB()
{
}

ReducedSteadyBB::ReducedSteadyBB(SteadyBB& FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC_t    = problem->inletIndexT.rows();
    N_BC      = problem->inletIndex.rows();
    Nphi_u    = problem->B_matrix.rows();
    Nphi_prgh = problem->K_matrix.cols();
    Nphi_t    = problem->Y_matrix.rows();
Info << "Bug1" << endl;
    // Create locally the velocity modes
    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        LUmodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        LUmodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        LUmodes.append(problem->supmodes[k]);
    }

    // Create locally the temperature modes including BC with liftfield
    for (label k = 0; k < problem->liftfieldT.size(); k++)
    {
        LTmodes.append(problem->liftfieldT[k]);
    }

    for (label k = 0; k < problem->NTmodes; k++)
    {
        LTmodes.append(problem->Tmodes[k]);
    }

    newton_object_sup = newton_steadyBB_sup(Nphi_u + Nphi_prgh + Nphi_t,
                        Nphi_u + Nphi_prgh + Nphi_t,
                        FOMproblem);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //
//Operator to evaluate the residual for the supremizer approach
int newton_steadyBB_sup::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_prgh);
    Eigen::VectorXd cTmp(Nphi_t);
    aTmp = x.head(Nphi_u);
    bTmp = x.segment(Nphi_u, Nphi_prgh);
    cTmp = x.tail(Nphi_t);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Diffusive Term
    Eigen::VectorXd m1 = problem->B_matrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
    // Continuity
    Eigen::VectorXd m3 = problem->P_matrix * aTmp;
    // Buoyancy Term
    Eigen::VectorXd m10 = problem->H_matrix * cTmp;
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // diffusive term temperature
    Eigen::VectorXd m6 = problem->Y_matrix * cTmp * (nu / Pr) ;
    for (label i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp;
        fvec(i) = m1(i)  - cc(0, 0) - m10(i) - m2(i);
    }

    for (label j = 0; j < Nphi_prgh; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    for (label j = 0; j < Nphi_t; j++)
    {
        label k = j + Nphi_u + Nphi_prgh;
        qq = aTmp.transpose() * problem->Q_matrix[j]  * cTmp;
        fvec(k) = m6(j) - qq(0, 0);
    }

    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - bc(j);
    }

    for (label j = 0; j < N_BC_t; j++)
    {
        label k = j + Nphi_u + Nphi_prgh;
        fvec(k) = x(k) - bc_t(j);
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_steadyBB_sup::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_steadyBB_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
Eigen::MatrixXd ReducedSteadyBB::solveOnline_sup(Eigen::MatrixXd& temp,
        Eigen::MatrixXd& vel, int param)
{
    vel_now = vel;
    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[0], LUmodes);
    y.tail(Nphi_t) = ITHACAutilities::get_coeffs(problem->Tfield[0], LTmodes);

    // Set size of online solution
    online_solutiont.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);

    if (problem->bcMethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }

        for (label j = 0; j < N_BC_t; j++)
        {
            label k = j + Nphi_prgh + Nphi_u;
            y(k) = temp(j, 0);
        } 
    }

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_steadyBB_sup> hnls(newton_object_sup);

    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.Pr = Pr;
    newton_object_sup.bc_t.resize(N_BC_t);
    newton_object_sup.bc.resize(N_BC);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC_t; j++)
    {
        newton_object_sup.bc_t(j) = temp(j, 0);
    }

    for (label j = 0; j < N_BC; j++)
    {
        newton_object_sup.bc(j) = vel_now(j, 0);
    }
    
    hnls.solve(y);  
    Eigen::VectorXd res(y);
    newton_object_sup.operator()(y, res);

    if (Pstream::master())
    {
        std::cout << "Solving for the parameter: " << param << std::endl;
    }

    if (res.norm() < 1e-5 && Pstream::master())
    {
        std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    else if (Pstream::master())
    {
        std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }

   
    online_solutiont(0,0) = param + 1;
    online_solutiont.col(0).tail(y.rows()) = y;

    return online_solutiont;
}


void ReducedSteadyBB::reconstruct_sup(fileName folder, int printEvery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            volVectorField uRec("uRec", LUmodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                uRec += LUmodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(uRec, name(online_solution[i](0, 0)), folder);
            volScalarField TRec("TRec", LTmodes[0] * 0);

            for (label j = 0; j < Nphi_t; j++)
            {
                TRec += LTmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(TRec, name(online_solution[i](0, 0)), folder);
            nextWrite += printEvery;
            UREC.append(uRec);
            TREC.append(TRec);
        }

        counter++;
    }
}




// ************************************************************************* //

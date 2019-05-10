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
Class
    reducedSteadyNS
Description
    A reduced problem for the stationary NS equations
SourceFiles
    reducedSteadyNS.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedSteadyNS class.
/// \dir
/// Directory containing the header and source files for the reducedSteadyNS class

#ifndef reducedPisoUnsteadyNS_H
#define reducedPisoUnsteadyNS_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "reducedProblem.H"
#include "unsteadyNS_piso.H"
#include "reducedUnsteadyNS.H"
#include "ITHACAutilities.H"
#include "EigenFunctions.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "steadyNS.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "ITHACAforces.H"
#include "volFields.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "Modes.H"
#include <iostream>

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> steady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class reducedPisoUnsteadyNS: public reducedProblem
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedPisoUnsteadyNS();

        /// Construct Null
        ///
        /// @param[in]  problem  a full order steadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit reducedPisoUnsteadyNS(unsteadyNS_piso& problem);

        ~reducedPisoUnsteadyNS() {};

        // Functions

        /// Scalar to store the current time
        scalar time;

        /// Scalar to store the time increment
        scalar dt;

        /// Scalar to store the final time if the online simulation
        scalar finalTime;

        /// Scalar to store the initial time if the online simulation
        scalar tstart;
        scalar nextwrite;
        scalar printevery;

	int counter2;

       // Eigen::VectorXd uresidual;
       // Eigen::VectorXd presidual;

      //  scalar U_norm_res(1);
      //  scalar P_norm_res(1);
   
      //  int counter = 0;

        /// Reconstructed pressure field
        PtrList<volScalarField> PREC;

        /// Recontructed velocity field
        PtrList<volVectorField> UREC;


	Eigen::MatrixXd online_solution;

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void solveOnline_Piso();

        ///
        /// @brief      Gets the Umatrix for the online solver: it just assembles the pure matrix.
        ///
        /// @param[in]      U     Velocity vector field
        /// @param[in]      p     Pressure scalar field
        ///
        /// @return     Velocity linear system.
        ///
        fvVectorMatrix get_Umatrix_Online(volVectorField& U, volScalarField& p);

        ///
        /// @brief      Gets the Pmatrix for the online solver: it just assembls the pure matrix.
        ///
        /// @param[in]      U     Velocity vector field
        /// @param[in]      p     Pressure scalar field
        ///
        /// @return     Pressure linear system.
        ///
        fvScalarMatrix get_Pmatrix_Online(volVectorField& U, volScalarField& p);

        ///
        /// @brief      It checks if the number of imposed boundary conditions is correct
        /// and set the inlet velocity equal to the given one.
        ///
        /// @param[in]  vel   Inlet velocity
        ///
        void setOnlineVelocity(Eigen::MatrixXd vel);

        /// Method to reconstruct the solutions from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printevery  The printevery
        ///
        void reconstruct_sup(fileName folder = "./ITHACAOutput/online_rec",
                             int printevery = 1);

        // Variables

        /// Lifted velocity modes.
        Modes<vector> ULmodes;

        /// Full problem.
        unsteadyNS_piso* problem;

        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;

        /// Counter.
        int counter = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif







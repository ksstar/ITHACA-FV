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
    reducedProblem
Description
    A reduced problem for the unsteady NS equations
SourceFiles
    reducedUnsteadyNSExplicit.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedUnsteadyNSExplicit class.
/// \dir
/// Directory containing the header and source files for the reducedUnsteadyNSExplicit class

#ifndef ReducedUnsteadyNSExplicit_H
#define ReducedUnsteadyNSExplicit_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedSteadyNS.H"
#include "unsteadyNSExplicit.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class reducedUnsteadyNSExplicit: public reducedSteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedUnsteadyNSExplicit();
        /// Construct Null
        ///
        /// @param      problem  a full order unsteadyNSExplicit problem
        ///
        explicit reducedUnsteadyNSExplicit(unsteadyNSExplicit& problem);

        ~reducedUnsteadyNSExplicit() {};

        /// Scalar to store the current time
        scalar time;

        /// Scalar to store the time increment
        double dt;

	/// number of phi modes
	int Nphi_phi;

        /// Scalar to store the final time if the online simulation
        scalar finalTime;

        /// Scalar to store the initial time if the online simulation
        scalar tstart;

        /// A variable for storing the reduced coefficients
        double storeEvery;

        /// A variable for exporting the fields
        double exportEvery;

        /// A variable for starting solving the reduced system from zero or not
        bool startFromZero = false;

        /// Pointer to the FOM problem
        unsteadyNSExplicit* problem;

        // Functions

        /// Method to perform an online solve using the EXPLICIT method
        ///
        /// @param[in]  vel_now   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snapshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline(Eigen::MatrixXd vel_now, label startSnap = 0);


        /// Method to reconstruct a solution from an online solve.
        ///
        /// @param[in]  folder      The folder where to output the solutions
        ///
        void reconstruct(fileName folder = "./online_rec");


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif


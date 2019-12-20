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
    ReducedSteadyBBTurb
Description
    A reduced problem for the stationary turbulent NS equations
SourceFiles
    ReducedSteadyBBTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedSteadyBBTurb class.
/// \dir
/// Directory containing the header and source files for the ReducedSteadyBBTurb class

#ifndef ReducedSteadyBBTurb_H
#define ReducedSteadyBBTurb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "ReducedSteadyNS.H"
#include "ReducedSteadyBB.H"
#include "SteadyBBTurb.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


/// Structure to implement a newton object for a stationary NS problem
struct newtonSteadyBBTurb: public newton_argument<double>
{
    public:
        newtonSteadyBBTurb() {}
        newtonSteadyBBTurb(int Nx, int Ny,
                           SteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            nphiNut(problem.nNutModes),
	    nphiTke(problem.nTkeModes),
            Nphi_p(problem.NPmodes),
	    Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
            N_BC_t(problem.inletIndexT.rows()),
            N_BC(problem.inletIndex.rows()),
            Nphi_prgh(problem.NPrghmodes),
            gNut(problem.nNutModes),
	    kNut(problem.nTkeModes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        SteadyBBTurb* problem;
        int Nphi_u;
        int nphiNut;
	int nphiTke;
        int Nphi_p;
        int N_BC;
	int Nphi_t;
        int N_BC_t;
        int Nphi_prgh;
        Eigen::VectorXd gNut;
	Eigen::VectorXd kNut;
        scalar nu;
        scalar Pr;
	scalar Prt;
 	Eigen::VectorXd bc_t;
        Eigen::VectorXd bc;
	Eigen::VectorXd pCnst;
	Eigen::VectorXd tauT;
	Eigen::VectorXd tauU;
	Eigen::VectorXd tauHF;
        Eigen::MatrixXd gradient_t;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/


/// Class where it is implemented a reduced problem for the <b> steady turbulent Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary turbulent Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedSteadyBBTurb: public ReducedSteadyBB
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedSteadyBBTurb();
        /// Construct Null
        ///
        /// @param      problem  a full order SteadyNSTurb problem
        ///
        ReducedSteadyBBTurb(SteadyBBTurb& problem);

        ~ReducedSteadyBBTurb() {};

	/// The matrix of the eddy viscosity RBF interoplated coefficients
        List <Eigen::MatrixXd> rbfCoeffMat;

	/// List of POD modes for pressure
	PtrList<volScalarField> LPmodes;
	/// List of POD modes for shifted pressure
	PtrList<volScalarField> LPrghmodes;

        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

	/// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> kFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> kModes;

        /// Reconstructed eddy viscosity field
        PtrList<volScalarField> nutRec;

        /// Vector of eddy viscosity RBF interoplated coefficients
        Eigen::VectorXd rbfCoeff;

	/// Vector of turbulent kinetic energy RBF interoplated coefficients
	Eigen::VectorXd tkeRbfCoeff;

	/// To determine index
	Eigen::VectorXd par_on_BC;

        /// Number of viscosity modes
        int nphiNut;

	/// Number of kinetic energy modes
	int nphiTke;

	/// Turbulent Prandtl number
 	scalar Prt;

        /// Newton Object to solve the nonlinear problem
        newtonSteadyBBTurb newtonObject;

        /// Pointer to the FOM problem
        SteadyBBTurb* problem;

        ///
        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnlineSUP(Eigen::MatrixXd velNow, Eigen::MatrixXd temp, int param = 0);

        /// Method to reconstruct the solutions from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printEvery  The printEvery
        ///
        void reconstructSUP(fileName folder = "./ITHACAOutput/online_rec",
                            int printEvery = 1);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif







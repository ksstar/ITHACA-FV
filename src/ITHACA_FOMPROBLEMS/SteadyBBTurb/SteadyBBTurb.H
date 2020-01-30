/*---------------------------------------------------------------------------*\
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
Class
    SteadyBBTurb
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    SteadyBBTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the SteadyBB class.
/// \dir
/// Directory containing the header and source files for the SteadyBBTurb class.

#ifndef SteadyBBTurb_H
#define SteadyBBTurb_H
#include "SteadyBB.H"
#include "SteadyNSTurb.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyBBTurb Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady Boussinesq problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a steady Boussinesq problem (two way coupling between momentum and energy equation)
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class SteadyBBTurb: public SteadyBB
{
    public:
        // Constructors
        /// Null constructor
        SteadyBBTurb();

        /// Construct with argc and argv
        SteadyBBTurb(int argc, char* argv[]);
        //~SteadyBBTurb() {};

        // Member Functions
	/// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;
	
	/// List of snapshots for the solution for thermal diffusivity
        PtrList<volScalarField> alphatField;

	/// List of snapshots for the solution for kinetic energy
	PtrList<volScalarField> kFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

	/// List of POD modes for thermal diffusivity
        PtrList<volScalarField> alphatModes;

	/// List of POD modes for kinetic energy
	PtrList<volScalarField> kModes;

        /// Create a samples for interpolation
        std::vector<SPLINTER::DataTable*> samples;

        /// Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplines;

        /// Turbulent viscosity matrix
        Eigen::MatrixXd btMatrix;

	 /// Turbulent viscosity matrix
        Eigen::MatrixXd K2_matrix;

	 /// Turbulent viscosity matrix
        Eigen::MatrixXd H_matrix4;
	Eigen::MatrixXd H_matrix5;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct2Matrix;
        Eigen::Tensor<double, 3 > ct2Tensor;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct1Matrix;
        Eigen::Tensor<double, 3 > ct1Tensor;

        /// Total Turbulent tensor
        List <Eigen::MatrixXd> cTotalMatrix;
        Eigen::Tensor<double, 3 > cTotalTensor;

        /// Total B Matrix
        Eigen::MatrixXd bTotalMatrix;

	/// Turbulent diffusivity term
        List <Eigen::MatrixXd> S_matrix;

	/// Turbulent diffusivity term
        List <Eigen::MatrixXd> S2_matrix;


        /// The matrix of L2 projection coefficients for the eddy viscosity
        Eigen::MatrixXd coeffL2;
	
	/// The matrix of L2 projection coefficients for the kinetic energy
	Eigen::MatrixXd tkeCoeffL2;

	/// The matrix of L2 projection coefficients for the kinetic energy
	Eigen::MatrixXd alphatCoeffL2;

        /// The vector of L2 projection coefficients for the eddy viscosity snapshot
        Eigen::VectorXd nutCoeff;

	/// The vector of L2 projection coefficients for the kinetic energy snapshot
	Eigen::VectorXd kCoeff;

        /// The way to compute the eddy viscosity snapshots
        word viscCoeff;

        /// Number of viscosity modes used for the projection
        label nNutModes;
	
	/// Number of kinetic energy modes used for the projection
	label nTkeModes;

	/// Number of kinetic energy modes used for the projection
	label nAlphatModes;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;
        autoPtr<volScalarField> _nuTilda;

	/// RBF functions radius
	double e = 1;

 	/// Kinetic energy field
	autoPtr<volScalarField> _k;

 	// Create a samples for interpolation
        std::vector<SPLINTER::DataTable*> tkeSamples;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> tkeRbfSplines;

	// Create a samples for interpolation
        std::vector<SPLINTER::DataTable*> alphatSamples;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> alphatRbfSplines;
        
	//--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now, fileName folder = "./ITHACAoutput/Offline");

        // Wrapped Proj. Methods;

        //--------------------------------------------------------------------------
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder used to save the reduced matrices.
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NPrgh, label NTmodes,
                        label NSUPmodes, label nNutModes, label nTkeModes, Eigen::MatrixXd para2);

	//--------------------------------------------------------------------------
        ///
        /// @brief      bt added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     bt matrix for turbulence treatment
        ///
        Eigen::MatrixXd btTurbulence(label NUmodes, label NSUPmodes);


        //--------------------------------------------------------------------------
        /// @brief      ct1 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm1(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct1 tensor for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor1(label NUmodes,
                label NSUPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm2(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor2(label NUmodes,
                label NSUPmodes, label nNutModes);

        ///
        /// @brief      S added matrix for the thermal turbulence treatement
        ///
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     S matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > temperature_turbulence_term(label NTmodes,
                label nNutModes);

	///
        /// @brief      S2 added matrix for the thermal turbulence treatement
        ///
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     S2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > temperature_turbulence_term2(label NTmodes,
                label nNutModes);

	//--------------------------------------------------------------------------
        /// Gradient of pressure
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the Gradient of pressure term.
        ///
        Eigen::MatrixXd pressure_gradient_term2(label NUmodes, label NPmodes,
                                               label NSUPmodes);

	Eigen::MatrixXd buoyant_term4(label NUmodes, label NTmodes, label NSUPmodes);
	Eigen::MatrixXd buoyant_term5(label NUmodes, label NTmodes, label NSUPmodes);

};

#endif











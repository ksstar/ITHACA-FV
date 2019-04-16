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
    unsteadyNS_simple
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    unsteadyNS_simple.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNS class.
/// \dir
/// Directory containing the header and source files for the unsteadyNS class.

#ifndef unsteadyNS_piso_H
#define unsteadyNS_piso_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "unsteadyNS_PISO.H"
#include "steadyNS.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "ITHACAforces.H"
#include "volFields.H"
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class unsteadyNS_piso Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class unsteadyNS_piso: public unsteadyNS_PISO
{

    public:
        // Constructors
        /// Null constructor
        unsteadyNS_piso();

        /// Construct with argc and argv
        unsteadyNS_piso(int argc, char* argv[]);
        ~unsteadyNS_piso() {};

	/// pisoControl
        //autoPtr<pisoControl> _piso;

	///dimensionedScalar nu;
        //autoPtr<dimensionedScalar> _nu;

        // Variables
        /// Initialization for the full velocity linear system.
        fvVectorMatrix* Ueqn_global;
        /// Initialization for the full pressure linear system.
        fvScalarMatrix* Peqn_global;

        //Functions
        ///
        /// @brief      Gets the offline velocity matrix and solve the system. Procedure of the PISO algorithm implemented
        /// in openFOAM is followed with same corrections and iterations.
        ///
        /// @param      U     Velocity field at previous step.
        /// @param      p     Pressure field at previous step.
        ///
        /// @return     The velocity matrix.
        ///
        fvVectorMatrix get_Umatrix(volVectorField& U, volScalarField& p);
        ///
        /// @brief      Gets the offline pressure matrix, solve the system and updates the residual for pressure.
        /// Procedure of the PISO algorithm implemented in openFOAM is followed with same corrections and iterations.
        ///
        /// @param      U     Velocity field at previous step.
        /// @param      p     Pressure field at previous step.
        ///
        /// @return     The velocity matrix.
        ///
        fvScalarMatrix get_Pmatrix(volVectorField& U, volScalarField& p);
        ///
        /// @brief      Offline solver for the whole Navier-Stokes problem.
        ///
        /// @param[in]  mu_now  Viscosity (parametrized) we want to solve the problem with.
        ///
        void truthSolve2(List<scalar> mu_now);
};

#endif











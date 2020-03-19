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
    unsteadyNS
Description
    Reduction class for a NON-Stationary NS plus energy equation problem
SourceFiles
    unsteadyNSExplicit.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNSExplicit class.
/// \dir
/// Directory containing the header and source files for the unsteadyNS class.

#ifndef unsteadyNSExplicit_H
#define unsteadyNSExplicit_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class unsteadyNSPISO Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem weakly coupled with the energy equation </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NST problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the unsteadyNS class */
class unsteadyNSExplicit: public steadyNS
{

    public:
        /// Constructors

        /// Construct Null
        unsteadyNSExplicit();

        /// Construct with argc and argv
        unsteadyNSExplicit(int argc, char* argv[]);
        //~unsteadyNST() {};

        // Specific variable for the unstationary case
        /// Start Time (initial time to start storing the snapshots)
        scalar startTime;

        /// Final time (final time of the simulation and consequently of the acquisition of the snapshots)
        scalar finalTime;

        /// Time step of the simulation
        scalar timeStep;

        /// Time step of the writing procedure
        scalar writeEvery = timeStep;

        /// Auxiliary variable to store the next writing instant
        scalar nextWrite;

        /// adjustTimeStep
        bool adjustTimeStep;

        /// maxCourant
        scalar maxCo;

        /// maxDeltaT
        scalar maxDeltaT;

        ///dimensionedScalar nu;
        autoPtr<dimensionedScalar> _nu;

	///dimensionedScalar dt;
        autoPtr<dimensionedScalar> _dt;

        // Counter used for the time-dependent boundary conditions
        label counter2 = 1;

        /// Time-dependent Boundary Method
        word timedepbcMethod;

        // Matrix containing the time-dependent boundary conditions for each time step
        Eigen::MatrixXd timeBCoff;

        // Time derivative numerical scheme order
        word timeDerivativeSchemeOrder;

	// Matrix to store norm of div(U) and div(phi)
	Eigen::MatrixXd DivNorm;

        /// pisoControl
        autoPtr<pisoControl> _piso;

        /// Reference pressure cell
        label pRefCell;

        /// Reference pressure value
        scalar pRefValue;

	scalar a21;
	scalar a31;
	scalar a32;
	scalar b1;
	scalar b2;
	scalar b3;
	double c1;
	double c2;
	double c3;

	/*volVectorfield U2;
	volVectorfield U3;
	volVectorfield U4;

	volVectorfield dU;


	volScalarfield p2;
	volScalarfield p3;
	volScalarfield p4;

	surfaceScalarField phi2;
	surfaceScalarField phi3;
	surfaceScalarField phi4;*/


        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// Time
        autoPtr<Time> _runTime;


        /// Perform a truthsolve
        void truthSolve(List<scalar> mu_now,
                        fileName folder = "./ITHACAoutput/Offline/");

        //--------------------------------------------------------------------------
        /// Function to check if the solution must be exported.
        ///
        /// @param      timeObject  The time object of OpenFOAM.
        ///
        /// @return     1 if we must write 0 elsewhere.
        ///
        bool checkWrite(Time& timeObject);


        /// Specific variable for the unstationary case

        
};

#endif




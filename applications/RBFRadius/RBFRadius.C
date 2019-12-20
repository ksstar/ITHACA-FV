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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    lift_and_drag

Author
    Giovanni Stabile, SISSA MathLab (International School for Advanced Studies) gstabile@sissa.it

Description
    Application to recover the lift and the drag after the simulation is performed

\*---------------------------------------------------------------------------*/

/// \file
/// \brief Application to recover the lift and the drag after the simulation is performed
/// \details In order to use this file one needs to prepare a FORCESdict file, in order to 
/// check the syntax one needs to check the \ref FORCESdict file.

/// \file FORCESdict
/// \brief Example of a FORCESdict file

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "ITHACAstream.H"
#include "ITHACAutilities.H"

#include <Eigen/Dense>

#include "Foam2Eigen.H"
#include <chrono>

#include "EigenFunctions.H"
#include <chrono>
#include <Eigen/SVD>
#include <Eigen/SparseLU>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

// #include "setRootCase.H"
// #include "createTime.H"
// #include "createMesh.H"

	std::string app = argv[1];
        word toRun("./radii"); // the offline samples.
        Eigen::MatrixXd radiiToRun = ITHACAstream::readMatrix(toRun);
        int N_run = radiiToRun.rows();

        for (label k = 0; k < N_run; k++)
        {

	   word removePOD = "rm -r ITHACAoutput/POD/";
	    word removeMatrices = "rm -r ITHACAoutput/Matrices/";
            int status6 = system(removePOD);
	    int status7 = system(removeMatrices);
            word radiusCommand = "sed -i 's/^RBFradius .*$/RBFradius " + name(radiiToRun(k,0)) + ";/' ./system/ITHACAdict";
            int status0 = system(radiusCommand);


            int status1 = system(app);

            word removeWeights = "rm -r ITHACAoutput/weightsSUP/";
            word fileCommand1 = "cp -r ./ITHACAoutput/rbfCoeffMat ./ITHACAoutput/CrossValidation"
            + name(k+1);

            int status2 = system(removeWeights);
            int status5 = system(fileCommand1);

        }
        return 0;
    }


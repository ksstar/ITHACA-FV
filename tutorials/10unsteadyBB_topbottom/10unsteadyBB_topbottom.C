/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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

Description
    Example of Boussinesq approximation for two way coupling NS-momentum equation
    and heat transport equation for enclosed flows.
SourceFiles
    10unsteadyBB_enclosed.C
\*---------------------------------------------------------------------------*/

#include "unsteadyBB.H"
#include "ITHACAPOD.H"
#include "reducedUnsteadyBB.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorial10top: public unsteadyBB
{
    public:
        explicit tutorial10top(int argc, char* argv[])
            :
            unsteadyBB(argc, argv),
            U(_U()),
            p(_p()),
            p_rgh(_p_rgh()),
            T(_T())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& p_rgh;
        volScalarField& T;

        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (label k = 0; k < par_BC.rows(); k++)
                {
                    for (label j = 0; j < par_BC.cols(); j++)
                    {
                        for (label i = 0; i < mu.cols(); i++)
                        {
                            mu_now[0] = mu(0, i);
                        }

                        assignBC(T, inletIndexT(j, 0), par_BC(k, j));
                    }

                    truthSolve(mu_now);
                }
            }
        }


        void onlineSolveFull(Eigen::MatrixXd par_BC, label para_set_BC, fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
            }
            else
            {
                mkDir(folder);
                ITHACAutilities::createSymLink(folder);
                label i = para_set_BC;

                for (label j = 0; j < par_BC.cols(); j++)
                {
                    assignBC(T, inletIndexT(j, 0), par_BC(i, j));
                }

                truthSolve(folder);
            }
        }

        void onlineSolveRead(fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Tfield_on, T, folder);
            }
            else
            {
            }
        }


        // Method to compute the lifting function for temperature
        void liftSolveT()
        {
            for (label k = 0; k < inletIndexT.rows(); k++)
            {
                Time& runTime = _runTime();
                fvMesh& mesh = _mesh();
                volScalarField T = _T();
                simpleControl simple(mesh);
                volScalarField& alphat = _alphat();
                dimensionedScalar& Pr = _Pr();
                dimensionedScalar& Prt = _Prt();
                label BCind = inletIndexT(k, 0);
                volScalarField Tlift("Tlift" + name(k), T);
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                Info << "Solving a lifting Problem" << endl;
                scalar t1 = 1;
                scalar t0 = 0;

                for (label j = 0; j < T.boundaryField().size(); j++)
                {
                    assignIF(Tlift, t0);

                    if (j == BCind)
                    {
                        assignBC(Tlift, j, t1);
                    }
                    else if (T.boundaryField()[BCind].type() == "fixedValue")
                    {
                        assignBC(Tlift, j, t0);
                    }
                    else
                    {
                    }
                }

                while (simple.correctNonOrthogonal())
                {
                    alphat = turbulence->nut() / Prt;
                    alphat.correctBoundaryConditions();
                    volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(alphaEff, Tlift)
                    );
                    TEqn.solve();
                    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                         << nl << endl;
                }

                Tlift.write();
                liftfieldT.append(Tlift);
            }
        }
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial10top example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./par_offline_BC");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online_BC("./par_online_BC");
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesTproj   = para.ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    /// Set the parameter ranges
    example.mu_range(0, 0) = 0.00001;
    example.mu_range(0, 1) = 0.00001;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    /// Set the inlet Temperature boundaries where there are non homogeneous boundary conditions
    example.inletIndexT.resize(2, 1);
    example.inletIndexT << 0, 1;
    /// Time parameters
    example.startTime = 0.0;
    example.finalTime = 10.0;
    example.timeStep = 0.005;
    example.writeEvery = 0.01;
    // Perform the Offline Solve;
    example.offlineSolve(par_off_BC);
    // Search the lift function for the temperature
    example.liftSolveT();
    // Create homogeneous basis functions for temperature
    example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
    // Perform a POD decomposition for velocity temperature and pressure fields
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example.podex, 0, 0,
                        NmodesOut);
    // Get reduced matrices
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesTproj,
                       NmodesSUPproj);
    // Resize the modes for projection
    example.Tmodes.resize(NmodesTproj);
    example.Umodes.resize(NmodesUproj);
    // Online part
    reducedUnsteadyBB reduced(example);
    // Set values of the online solve
    reduced.nu = 0.00001;
    reduced.Pr = 0.71;
    reduced.tstart = 0.0;
    reduced.finalTime = 10;
    reduced.dt = 0.005;
    // No parametrization of velocity on boundary
    Eigen::MatrixXd vel_now_BC(0, 0);

    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (par_on_BC.rows()); k++)
    {
        Eigen::MatrixXd temp_now_BC(2, 1);
        temp_now_BC(0, 0) = par_on_BC(k, 0);
        temp_now_BC(1, 0) = par_on_BC(k, 1);
        reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 2);
    }
    // for which the offline solve has been performed
    example.onlineSolveRead("./ITHACAoutput/Offline/");
    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::error_listfields(
                                         example.Tfield_on, reduced.TREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");
    exit(0);
}

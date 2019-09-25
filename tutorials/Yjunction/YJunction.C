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
Description
    Example of an unsteady NS Reduction Problem
SourceFiles
    YJunction.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "reducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorialY: public unsteadyNS
{
    public:
        explicit tutorialY(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
		    U = U0;
		    p = P0;
                    mu_now[0] = mu(0, 0);
                    truthSolve(mu_now);
	            restart();
            }
        }

	void onlineSolveFull( fileName folder)
        {
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (ITHACAutilities::check_folder(folder))
            {
            }
            else
            {
                mkDir(folder);
                ITHACAutilities::createSymLink(folder);

		U = U0;
		p = P0;
                mu_now[0] = mu(0, 0);
                truthSolve(mu_now, folder);
	        restart();
            }
        }

        void onlineSolveRead(fileName folder)
        {
            List<scalar> mu_now(1);
	    volVectorField U0 = U;
	    volScalarField P0 = p;

            if (ITHACAutilities::check_folder(folder))
            {
		ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Pfield_on, p, folder);
            }
            else
	    {
            }
        }

// Method to compute the lifting functions for this tutorial
void liftSolve()
{
    for (label k = 0; k < 1; k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[0] = 1; // velocity magnitude = 1  and inlet 1
        v1[1] = -1; // velocity magnitude = 1  and inlet 1
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(Ulift);
    }



    for (label k = 1; k < 2; k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[0] = 1; // velocity magnitude = 1  and inlet 2
        v1[1] = 1; // velocity magnitude = 1  and inlet 2
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(Ulift);
    }


}

};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorialY object
    tutorialY example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./timeBCoffline");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    // the samples which will be used for setting the boundary condition in the online stage
   // word par_online_BC("./timeBConline");
   // Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read parameters from ITHACAdict file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    int NmodesOut = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 20);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();

    if (example.bcMethod == "lift")
    {
        example.inletIndex.resize(2, 2); 
        example.inletIndex(0, 0) = 1;  // Patch inlet 1
        example.inletIndex(0, 1) = 0;  
        example.inletIndex(1, 0) = 2;  // Patch inlet 2
        example.inletIndex(1, 1) = 0;  
   }
   else if (example.bcMethod == "penalty")
   {
        // Set the inlet boundaries where we have non homogeneous boundary conditions
        example.inletIndex.resize(4, 2); // rows: total number of patches 
        example.inletIndex(0, 0) = 1;  // Patch inlet 1
        example.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
        example.inletIndex(1, 0) = 1;  // Patch inlet 1: y-direction
        example.inletIndex(1, 1) = 1;  // Patch inlet 2
        example.inletIndex(2, 0) = 2;  // Patch inlet 2: x-direction
        example.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
        example.inletIndex(3, 0) = 2;  // Patch inlet 2: x-direction
        example.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
    }

    // Time parameters
    example.startTime = 0;
    example.finalTime = 12;
    example.timeStep = 0.0005;
    example.writeEvery = 0.03;
    example.Dim = 2;

    int totalPeriods = 4;

    /* // Period 1: for 0 IC to steady state for U1=V1=U2=V2 = 1 m/s --> 3seconds
    Eigen::VectorXd option1 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.70711,-0.70711);
    Eigen::VectorXd option1b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.70711,0.70711);
    Eigen::VectorXd option1c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.35355,-0.35335);
    Eigen::VectorXd option1d = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.35355,0.35355);
    // Period 2: Linear decrease to 50% --> 3seconds
    Eigen::VectorXd option2 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.70711,-0.35355);
    Eigen::VectorXd option2b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.70711,0.35355);
    Eigen::VectorXd option2c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.70711,0.5);
    Eigen::VectorXd option2d = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.70711,0.5);
    // Period 3: Linear increase to 50% --> 3seconds);
    Eigen::VectorXd option3 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.5,0.70711);
    Eigen::VectorXd option3b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.35355,-0.70711);
    Eigen::VectorXd option3c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.35355,0.70711);


    example.timeBCoff.resize(example.inletIndex.rows(), option1.size()*totalPeriods);

if (example.bcMethod == "penalty" && example.timedepbcMethod == "yes")
    {   
         //1st Period
         example.timeBCoff.row(0).head(option1.size()) = option1b.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).head(option1.size()) = option1.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).head(option1.size()) = option2b.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).head(option1.size()) = option2b.col(0);  //Patch inlet 4: x-direction

         //2nd Period
         example.timeBCoff.row(0).segment(option1.size(),option1.size()) = option2b.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).segment(option1.size(),option1.size()) = option2.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).segment(option1.size(),option1.size()) = option1d.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).segment(option1.size(),option1.size()) = option1d.col(0);  //Patch inlet 4: x-direction

         //3rd Period
         example.timeBCoff.row(0).segment(2*option1.size(),option1.size()) = option1d.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).segment(2*option1.size(),option1.size()) = option1c.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).segment(2*option1.size(),option1.size()) = option3c.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).segment(2*option1.size(),option1.size()) = option3c.col(0);  //Patch inlet 4: x-direction

         //4th Period
         example.timeBCoff.row(0).tail(option1.size()) = option3c.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).tail(option1.size()) = option3b.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).tail(option1.size()) = option1b.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).tail(option1.size()) = option1b.col(0);  //Patch inlet 4: x-direction

std::cout << "################## Here ##################" << std::endl;

ITHACAstream::exportMatrix(example.timeBCoff, "timeBCoff", "eigen",
                               "./ITHACAoutput/timeBCoff");

}
else
{
} */

Eigen::MatrixXd vel_now;
    vel_now = example.timeBCoff;

    // Parameter set u2A:
    Eigen::VectorXd VEL2Ax = Eigen::VectorXd::LinSpaced(3.0/example.timeStep+1,-0.707110,-0.473760);
    Eigen::VectorXd VEL2Ay = Eigen::VectorXd::LinSpaced(3.0/example.timeStep+1,0.707110,0.473760); // 1 to 0.67
    Eigen::VectorXd VEL2AL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep+1,0.707110,0.670000);

    // Parameter set u2B:
    Eigen::VectorXd VEL2Bx = Eigen::VectorXd::LinSpaced(4.0/example.timeStep,-0.47376,-0.70711);
    Eigen::VectorXd VEL2By = Eigen::VectorXd::LinSpaced(4.0/example.timeStep,0.47376,0.70711);
    Eigen::VectorXd VEL2BL = Eigen::VectorXd::LinSpaced(4.0/example.timeStep,0.67,0.70711);     // 0.67 to 1

    // Parameter set u2C:
    Eigen::VectorXd VEL2Cx = Eigen::VectorXd::LinSpaced(5.0/example.timeStep,-0.707110,-0.353550);
    Eigen::VectorXd VEL2Cy = Eigen::VectorXd::LinSpaced(5.0/example.timeStep,0.707110,0.353550);
    Eigen::VectorXd VEL2CL = Eigen::VectorXd::LinSpaced(5.0/example.timeStep,0.707110,0.50);   // 1 to 0.5

    // Parameter set u2D:
    Eigen::VectorXd VEL2Dx = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,-0.35355,-0.42426);
    Eigen::VectorXd VEL2Dy = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.35355,0.42426);
    Eigen::VectorXd VEL2DL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.5,0.6); // 0.5 to 0.6

    // Parameter set u2E:
    Eigen::VectorXd VEL2Ex = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,-0.42426,-0.42426);
    Eigen::VectorXd VEL2Ey = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,0.42426,0.42426);
    Eigen::VectorXd VEL2EL = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,0.6,0.6); // 0.6 to 0.6

    // Parameter set u2F:
    Eigen::VectorXd VEL2Fx = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,-0.42426,-0.50912);
    Eigen::VectorXd VEL2Fy = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.42426,0.50912);
    Eigen::VectorXd VEL2FL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.6,0.72); // 0.6 to 0.6

    // Parameter set u1A:
    Eigen::VectorXd VEL1Ax = Eigen::VectorXd::LinSpaced(2.0/example.timeStep+1,-0.707110,-0.671750);
    Eigen::VectorXd VEL1AL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep+1,0.707110,0.671750);
//    Eigen::VectorXd VEL1AL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep+1,1,0.95); // 1 to 0.95

    // Parameter set u1B:
    Eigen::VectorXd VEL1Bx = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,-0.67175,-0.50912);
    Eigen::VectorXd VEL1BL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.67175,0.50912);
   // Eigen::VectorXd VEL1BL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.95,0.72); // 0.95 to 0.72

    // Parameter set u1C:
    Eigen::VectorXd VEL1Cx = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,-0.50912,-0.50912);
    Eigen::VectorXd VEL1CL = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,0.50912,0.50912); 
  //  Eigen::VectorXd VEL1CL = Eigen::VectorXd::LinSpaced(1.0/example.timeStep,0.72,0.72);

    // Parameter set u1D:
    Eigen::VectorXd VEL1Dx = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,-0.50912,-0.35355);
    Eigen::VectorXd VEL1DL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.50912,0.35355);
  //  Eigen::VectorXd VEL1DL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.72,0.5);

    // Parameter set u1E:
    Eigen::VectorXd VEL1Ex = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,-0.35355,-0.58690);
    Eigen::VectorXd VEL1EL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.35355,0.58690);
  //  Eigen::VectorXd VEL1EL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.5,0.83);

    // Parameter set u1F:
    Eigen::VectorXd VEL1Fx = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,-0.58690,-0.63640);
    Eigen::VectorXd VEL1FL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.58690,0.63640);
 //   Eigen::VectorXd VEL1FL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.83,0.90);

    // Parameter set u1G:
    Eigen::VectorXd VEL1Gx = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,-0.63640,-0.50912);
    Eigen::VectorXd VEL1GL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.63640,0.50912);
 //   Eigen::VectorXd VEL1GL = Eigen::VectorXd::LinSpaced(2.0/example.timeStep,0.90,0.72);

// Parameter set u1G:
    Eigen::VectorXd VEL1Hx = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,-0.50912,-0.67175);
    Eigen::VectorXd VEL1HL = Eigen::VectorXd::LinSpaced(3.0/example.timeStep,0.50912,0.67175);



vel_now.resize(example.inletIndex.rows(), 36001);

    if (example.bcMethod == "penalty" && example.timedepbcMethod == "yes")
    {   
        

        // section A
    	vel_now.row(0).head(VEL1Ax.size()) = VEL1AL.col(0);  
    	vel_now.row(1).head(VEL1Ax.size()) = VEL1Ax.col(0);  
    	vel_now.row(2).head(VEL2Ax.size()) = VEL2Ay.col(0); 
    	vel_now.row(3).head(VEL2Ay.size()) = VEL2Ay.col(0);  
        // section B
        vel_now.row(0).segment(VEL1Ax.size(),VEL1Bx.size()) = VEL1BL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size(),VEL1Bx.size()) = VEL1Bx.col(0);  
        vel_now.row(2).segment(VEL2Ax.size(),VEL2Bx.size()) = VEL2By.col(0); 
        vel_now.row(3).segment(VEL2Ay.size(),VEL2By.size()) = VEL2By.col(0);  
        // section C
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size(),VEL1Cx.size()) = VEL1CL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size()+VEL1Bx.size(),VEL1Cx.size()) = VEL1Cx.col(0);  
        vel_now.row(2).segment(VEL2Ax.size()+VEL2Bx.size(),VEL2Cx.size()) = VEL2Cy.col(0); 
        vel_now.row(3).segment(VEL2Ay.size()+VEL2By.size(),VEL2Cy.size()) = VEL2Cy.col(0);  
        // section D
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size(),VEL1Dx.size()) = VEL1DL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size(),VEL1Dx.size()) = VEL1Dx.col(0);  
        vel_now.row(2).segment(VEL2Ax.size()+VEL2Bx.size()+VEL2Cx.size(),VEL2Dx.size()) = VEL2Dy.col(0); 
        vel_now.row(3).segment(VEL2Ay.size()+VEL2By.size()+VEL2Cy.size(),VEL2Dy.size()) = VEL2Dy.col(0);  
   	// section E
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size(),VEL1Ex.size()) = VEL1EL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size(),VEL1Ex.size()) = VEL1Ex.col(0);  
        vel_now.row(2).segment(VEL2Ax.size()+VEL2Bx.size()+VEL2Cx.size()+VEL2Dx.size(),VEL2Ex.size()) = VEL2Ey.col(0); 
        vel_now.row(3).segment(VEL2Ax.size()+VEL2Bx.size()+VEL2Cx.size()+VEL2Dx.size(),VEL2Ey.size()) = VEL2Ey.col(0);  
	// section F
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size(),VEL1Fx.size()) = VEL1FL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size(),VEL1Fx.size()) = VEL1Fx.col(0);
        vel_now.row(2).tail(VEL2Fx.size()) = VEL2Fy.col(0); 
        vel_now.row(3).tail(VEL2Fy.size()) = VEL2Fy.col(0);
	// section G
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size()+VEL1Fx.size(),VEL1Gx.size()) = VEL1GL.col(0);  
        vel_now.row(1).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size()+VEL1Fx.size(),VEL1Gx.size()) = VEL1Gx.col(0);  
        // section H
        vel_now.row(0).tail(VEL1Hx.size()) = VEL1HL.col(0);  
        vel_now.row(1).tail(VEL1Hx.size()) = VEL1Hx.col(0); 

std::cout << "################## Online Penalty ##################" << std::endl;

ITHACAstream::exportMatrix(vel_now, "timeBConPenalty", "eigen",
                               "./ITHACAoutput/timeBCoff");

    }
    else if (example.bcMethod == "lift" && example.timedepbcMethod == "yes")
    {
        //vel_now.resize(example.inletIndex.rows(), 36001);
        // section A
    	vel_now.row(0).head(VEL1Ax.size()) = VEL1AL.col(0);  
    	vel_now.row(1).head(VEL2Ay.size()) = VEL2Ay.col(0); 
        // section B
        vel_now.row(0).segment(VEL1Ax.size(),VEL1Bx.size()) = VEL1BL.col(0);  
        vel_now.row(1).segment(VEL2Ay.size(),VEL2By.size()) = VEL2By.col(0); 
        // section C
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size(),VEL1Cx.size()) = VEL1CL.col(0);  
        vel_now.row(1).segment(VEL2Ay.size()+VEL2By.size(),VEL2Cy.size()) = VEL2Cy.col(0);  
        // section D
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size(),VEL1Dx.size()) = VEL1DL.col(0);  
        vel_now.row(1).segment(VEL2Ay.size()+VEL2By.size()+VEL2Cy.size(),VEL2Dy.size()) = VEL2Dy.col(0);    
   	// section E
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size(),VEL1Ex.size()) = VEL1EL.col(0);  
        vel_now.row(1).segment(VEL2Ay.size()+VEL2By.size()+VEL2Cy.size()+VEL2Dy.size(),VEL2Ey.size()) = VEL2Ey.col(0); 
	// section F
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size(),VEL1Fx.size()) = VEL1FL.col(0); 
        vel_now.row(1).tail(VEL2Fy.size()) = VEL2Fy.col(0); 
	// section G
        vel_now.row(0).segment(VEL1Ax.size()+VEL1Bx.size()+VEL1Cx.size()+VEL1Dx.size()+VEL1Ex.size()+VEL1Fx.size(),VEL1Gx.size()) = VEL1GL.col(0); 
        // section H
        vel_now.row(0).tail(VEL1Hx.size()) = VEL1HL.col(0);


std::cout << "################## Online Lift ##################" << std::endl;

ITHACAstream::exportMatrix(vel_now, "timeBConControl", "eigen",
                               "./ITHACAoutput/timeBCoff");
    } 


exit(0); 

example.timeBCoff = par_off_BC;



    // Perform The Offline Solve;
    example.offlineSolve();

    // Perform POD
    auto start_POD = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "lift")
    {

        if (example.aveMethod == "mean")
        { 
            // Search the lift function
   	    example.liftSolve();
    	    // Create homogeneous basis functions for velocity
    	    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);

            // Compute the average velocity and temperature fields
            ITHACAutilities::getAverage(example.Uomfield, example.Umean, example.Usub);
            ITHACAutilities::getAverage(example.Pfield, example.Pmean, example.Psub);

    	    // Perform a POD decomposition for velocity and pressure
   	    ITHACAPOD::getModes(example.Usub, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    	    ITHACAPOD::getModes(example.Psub, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);

        }
	else
	{
            // Search the lift function
   	    example.liftSolve();
  	    // Normalize the lifting function
	    ITHACAutilities::normalizeFields(example.liftfield);
    	    // Create homogeneous basis functions for velocity
    	    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    	    // Perform a POD decomposition for velocity and pressure
   	    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    	    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
        }
    }
    else if (example.bcMethod == "penalty")
    {
        if (example.aveMethod == "mean")
        { 
            // Compute the average velocity and temperature fields
            ITHACAutilities::getAverage(example.Ufield, example.Umean, example.Usub);
            ITHACAutilities::getAverage(example.Pfield, example.Pmean, example.Psub);

	    // Perform a POD decomposition for velocity and pressure
   	    ITHACAPOD::getModes(example.Usub, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    	    ITHACAPOD::getModes(example.Psub, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);

        }
	else
	{
            // Perform a POD decomposition for velocity and pressure
            ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
            ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
	}
    }

    auto finish_POD = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
    std::cout << "elapsed_POD: " << elapsed_POD.count() << " seconds.";
    std::cout << std::endl;


    // PPE method
    NmodesSUPproj = 0;

    // Solve the supremizer problem
/*    auto start_sup = std::chrono::high_resolution_clock::now();
    example.solvesupremizer("modes");
    auto finish_sup = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sup = finish_sup - start_sup;
    std::cout << "elapsed_sup: " << elapsed_sup.count() << " seconds.";
    std::cout << std::endl; */



    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl; 

/*
 // Create a list with number of modes for which the projection needs to be performed
    Eigen::MatrixXd List_of_modes(NmodesOut-0, 1);
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        List_of_modes(i, 0) = i + 1;
    }

    // Export with number of modes for which the projection needs to be performed
    ITHACAstream::exportMatrix(List_of_modes, "List_of_modes", "eigen", "./ITHACAoutput/l2error");

    PtrList<volVectorField> Umodes;

    if (example.liftfield.size() != 0)
    {
        for (label k = 0; k < example.liftfield.size(); k++)
        {
            Umodes.append(example.liftfield[k]);
        }
    }

    if (example.aveMethod == "mean")
    {
      Umodes.append(example.Umean[0]);
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
            Umodes.append(example.Umodes[k]);
    }

    if (example.supex == 1)
    {
        for (label k = 0; k < NmodesSUPproj; k++)
        {
            Umodes.append(example.supmodes[k]);
        }
    }

    PtrList<volScalarField> Pmodes;

    if (example.aveMethod == "mean")
    {
       Pmodes.append(example.Pmean[0]);
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        Pmodes.append(example.Pmodes[k]);
    }
   


    // Perform the projection for all number of modes in list List_of_modes
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size()-1, List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixP(example.Pfield.size()-1, List_of_modes.rows());

    // Calculate the coefficients and L2 error and store the error in a matrix for each number of modes
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
	std::cout << "Number of modes = " << i << std::endl;

	if (example.aveMethod == "mean")
        { 
	    Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Usub,Umodes,
                                   List_of_modes(i, 0) + example.liftfield.size() + example.Umean.size() + NmodesSUPproj);
            Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Psub, Pmodes,
                                 List_of_modes(i, 0) + example.Pmean.size());
            PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                		 Umodes, coeffU, List_of_modes(i, 0) + example.liftfield.size() + example.Umean.size() + NmodesSUPproj);
            PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                		 Pmodes, coeffP, List_of_modes(i, 0) + example.Pmean.size());

            L2errorProjMatrixU.col(i) = ITHACAutilities::error_listfields_min_IC(example.Usub, rec_fieldU);
            L2errorProjMatrixP.col(i) = ITHACAutilities::error_listfields_min_IC(example.Psub, rec_fieldP);

	}
	else
	{
            Eigen::MatrixXd coeffU = ITHACAutilities::get_coeffs(example.Ufield,Umodes,
                                   List_of_modes(i, 0) + example.liftfield.size() + NmodesSUPproj);
            Eigen::MatrixXd coeffP = ITHACAutilities::get_coeffs(example.Pfield, Pmodes,
                                 List_of_modes(i, 0));
            PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstruct_from_coeff(
                		 Umodes, coeffU, List_of_modes(i, 0));
            PtrList<volScalarField> rec_fieldP = ITHACAutilities::reconstruct_from_coeff(
                		 Pmodes, coeffP, List_of_modes(i, 0));
            L2errorProjMatrixU.col(i) = ITHACAutilities::error_listfields_min_IC(example.Ufield, rec_fieldU);
            L2errorProjMatrixP.col(i) = ITHACAutilities::error_listfields_min_IC(example.Pfield, rec_fieldP);
	}
    } 



    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixP, "L2errorProjMatrixP", "eigen",
                               "./ITHACAoutput/l2error"); */

    // Resize the modes for projection
    example.Umodes.resize(NmodesUproj);
    example.Pmodes.resize(NmodesPproj); 


    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.01;
    reduced.tstart = 0;
    reduced.finalTime = 18;
    reduced.dt = 0.0005;


    
    reduced.maxIter = 100;
    reduced.tolerance = 1e-5;
    reduced.timeSteps = 5;

    auto start_penalty = std::chrono::high_resolution_clock::now();
    if (example.bcMethod == "penalty")
    {
            // set initial quess for penalty factors
	    reduced.tauInit = Eigen::MatrixXd::Zero(4,1);
            reduced.tauInit <<  1e-6,1e-6,1e-6,1e-6; 
            reduced.tauU = reduced.penalty_PPE_time(vel_now, reduced.tauInit);
            //reduced.tauU = reduced.tauInit;
	
    }
    auto finish_penalty = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_penalty = finish_penalty - start_penalty;
    std::cout << "elapsed_penalty: " << elapsed_penalty.count() << " seconds.";
    std::cout << std::endl;  

    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (1); k++)
    {
	auto start_ROM = std::chrono::high_resolution_clock::now();
        reduced.solveOnline_PPE(vel_now, k);
std::cout << "HERE" << std::endl;
        auto finish_ROM = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
        std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
        std::cout << std::endl;
//exit(0);
        auto start_ROM_REC = std::chrono::high_resolution_clock::now();
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionPPE", 60);
        auto finish_ROM_REC = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM_REC = finish_ROM_REC - start_ROM_REC;
        std::cout << "elapsed_ROM_REC: " << elapsed_ROM_REC.count() << " seconds.";
        std::cout << std::endl;
    }

  
   // Performing full order simulation for second parameter set 
    tutorialY HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.01;
    HFonline2.mu_range(0, 1) = 0.01;
    HFonline2.genEquiPar();
    HFonline2.inletIndex.resize(4, 2); // rows: total number of patches 
    HFonline2.inletIndex(0, 0) = 1;  // Patch inlet 1
    HFonline2.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
    HFonline2.inletIndex(1, 0) = 1;  // Patch inlet 1: y-direction
    HFonline2.inletIndex(1, 1) = 1;  // Patch inlet 2
    HFonline2.inletIndex(2, 0) = 2;  // Patch inlet 2: x-direction
    HFonline2.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
    HFonline2.inletIndex(3, 0) = 2;  // Patch inlet 2: x-direction
    HFonline2.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 18;
    HFonline2.timeStep = 0.0005;
    HFonline2.writeEvery = 0.03;
    HFonline2.Dim = 2;

    // Set values BCs for each timeSetp:
    //HFonline2.timeBCoff.resize(HFonline2.inletIndex.rows(), option1.size());
    HFonline2.timeBCoff = vel_now;

    // Reconstruct the online solution
    HFonline2.onlineSolveFull("./ITHACAoutput/HFonline2");

    // Reading high-fidelity solutions for the parameter set
    // for which the offline solve has been performed (skipping IC)
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields(
                                         example.Pfield_on, reduced.PREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2error");


//Post-Process
    Eigen::MatrixXd PostP(example.Ufield_on.size(), 2); 

    for (label i = 0; i < example.Ufield_on.size(); i++)
    {
	PostP(i, 0) =  0.5*fvc::domainIntegrate(example.Ufield_on[i] & example.Ufield_on[i]).value();

	PostP(i, 1) = 0.5*fvc::domainIntegrate(reduced.UREC[i] & reduced.UREC[i]).value();
	
    }

    ITHACAstream::exportMatrix(PostP, "PostP", "eigen", "./ITHACAoutput/PostProcess");


exit(0); 
} 



  



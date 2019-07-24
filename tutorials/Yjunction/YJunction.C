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
#include<math.h>
#include<iomanip>

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

        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            Vector<double> inl(0, 0, 0);
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
                //for (label i = 0; i < par_BC.rows(); i++)
                //{
		    U = U0;
		    p = P0;
                    //inl[inletIndex(0,i+1)] = par_BC(i, 0);
                    mu_now[0] = mu(0, 0);
                    //assignBC(U, inletIndex(i,0), inl);
                    truthSolve(mu_now);
	            restart();
                //}
            }
        }



	void onlineSolveFull(Eigen::MatrixXd par_BC, label para_set_BC, fileName folder)
        {
 	    Vector<double> inl(0, 0, 0);
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
                label k = para_set_BC;

	        for (label i = 0; i < par_BC.rows(); i++)
                {
		    U = U0;
		    p = P0;
                    inl[inletIndex(0,i+1)] = par_BC(i, k);
                    mu_now[0] = mu(0, 0);
                    assignBC(U, inletIndex(i,0), inl);
                    truthSolve(mu_now, folder);
	            restart();
                }

            }
        }

        void onlineSolveRead(fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Pfield_on, p, folder);
            }
            else
            {
            }
        }



// Method to compute the lifting function
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
        v1[0] = -1;
        v1[1] = -1;
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
        v1[0] = -1;
        v1[1] = 1;
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
    word par_offline("./par_offline");
    Eigen::MatrixXd par_off = ITHACAstream::readMatrix(par_offline);
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online("./par_online");
    Eigen::MatrixXd par_on = ITHACAstream::readMatrix(par_online);
    // Read parameters from ITHACAdict file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    int NmodesOut     = para.ITHACAdict->lookupOrDefault<int>("NmodesOut", 20);
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

   // if (example.bcMethod == "lift")
  //  {
    ///    example.inletIndex.resize(2, 2); 
   //     example.inletIndex(0, 0) = 2;  // Patch inlet 1
   //     example.inletIndex(0, 1) = 0;  
   //     example.inletIndex(1, 0) = 3;  // Patch inlet 2
  //      example.inletIndex(1, 1) = 0;  

   // }
   // else
  //  {
        // Set the inlet boundaries where we have non homogeneous boundary conditions
        example.inletIndex.resize(4, 2); // rows: total number of patches 
        example.inletIndex(0, 0) = 2;  // Patch inlet 1
        example.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
        example.inletIndex(1, 0) = 2;  // Patch inlet 1: y-direction
        example.inletIndex(1, 1) = 1;  // Patch inlet 2
        example.inletIndex(2, 0) = 3;  // Patch inlet 2: x-direction
        example.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
        example.inletIndex(3, 0) = 3;  // Patch inlet 2: x-direction
        example.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
  //  }

    // Time parameters
    example.startTime = 0;
    example.finalTime = 12;
    example.timeStep = 0.0005;
    example.writeEvery = 0.03;

    int totalPeriods = 4;
example.Dim = 2;

    // Period 1: for 0 IC to steady state for U1=V1=U2=V2 = 1 m/s --> 3seconds
    Eigen::VectorXd option1 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-1,-1);
    Eigen::VectorXd option1b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,1,1);
    Eigen::VectorXd option1c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.35355,-0.35335);
    Eigen::VectorXd option1d = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.35355,0.35355);
    // Period 2: Linear decrease to 50% --> 3seconds
    Eigen::VectorXd option2 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-1,-0.35355);
    Eigen::VectorXd option2b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,1,0.35355);
    Eigen::VectorXd option2c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,1,0.5);
   
    // Period 3: Linear increase to 50% --> 3seconds);
    Eigen::VectorXd option3 = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.5,1);
    Eigen::VectorXd option3b = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,-0.35355,-1);
    Eigen::VectorXd option3c = Eigen::VectorXd::LinSpaced((example.finalTime/totalPeriods)/example.timeStep+1,0.35355,1);


    example.timeBCoff.resize(example.inletIndex.rows(), option1.size()*totalPeriods);
    //example.timeBCoff.row(0) = option1.col(0);  //Patch inlet 1: x-direction
   // example.timeBCoff.row(1) = option1.col(0);  // Patch inlet 1: y-direction

   // if (example.bcMethod == "penalty")
   // {
   //      example.timeBCoff.row(2) = option1.col(0); //Patch inlet 2: x-direction
   //      example.timeBCoff.row(3) = option1b.col(0);  //Patch inlet 4: x-direction
   // }
 

   // if (example.bcMethod == "penalty" && example.timedepbcMethod == "yes")
  //  {  //1st Period
         example.timeBCoff.row(0).head(option1.size()) = option1.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).head(option1.size()) = option1.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).head(option1.size()) = option2.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).head(option1.size()) = option2b.col(0);  //Patch inlet 4: x-direction

         //2nd Period
         example.timeBCoff.row(0).segment(option1.size(),option1.size()) = option2.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).segment(option1.size(),option1.size()) = option2.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).segment(option1.size(),option1.size()) = option1c.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).segment(option1.size(),option1.size()) = option1d.col(0);  //Patch inlet 4: x-direction

         //3rd Period
         example.timeBCoff.row(0).segment(2*option1.size(),option1.size()) = option1c.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).segment(2*option1.size(),option1.size()) = option1c.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).segment(2*option1.size(),option1.size()) = option3b.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).segment(2*option1.size(),option1.size()) = option3c.col(0);  //Patch inlet 4: x-direction

         //4th Period
         example.timeBCoff.row(0).tail(option1.size()) = option3b.col(0);  //Patch inlet 1: x-direction
         example.timeBCoff.row(1).tail(option1.size()) = option3b.col(0);  // Patch inlet 1: y-direction
         example.timeBCoff.row(2).tail(option1.size()) = option1.col(0); //Patch inlet 2: x-direction
         example.timeBCoff.row(3).tail(option1.size()) = option1b.col(0);  //Patch inlet 4: x-direction


   // }




    // Perform The Offline Solve;
    example.offlineSolve(par_off);

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
                        NmodesOut);
            ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesOut);
	}
    }

    auto finish_POD = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_POD = finish_POD - start_POD;
    std::cout << "elapsed_POD: " << elapsed_POD.count() << " seconds.";
    std::cout << std::endl;

    // PPE method
    NmodesSUPproj = 0;

    // Reduced Matrices
    auto start_matrix = std::chrono::high_resolution_clock::now();
    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    auto finish_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_matrix = finish_matrix - start_matrix;
    std::cout << "elapsed_matrix: " << elapsed_matrix.count() << " seconds.";
    std::cout << std::endl; 

/*
// Create a list with number of modes for which the projection needs to be performed
    Eigen::MatrixXd List_of_modes(NmodesOut-5, 1);
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
                               "./ITHACAoutput/l2error"); 

    // Resize the modes for projection
    example.Umodes.resize(NmodesUproj);
    example.Pmodes.resize(NmodesPproj);*/

    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.01;
    reduced.tstart = 0;
    reduced.finalTime = 12;
    reduced.dt = 0.0005;
    reduced.maxIter = 500;
    reduced.tolerance = 1e-7;
    reduced.timeSteps = 3;


    Eigen::MatrixXd vel_now;

    if (example.bcMethod == "lift")
    {
         vel_now.resize(2,1);
         vel_now(0,0) = 1;
         vel_now(1,0) = 1;
    }
    else
    {
         vel_now = example.timeBCoff;
    }

    if (example.bcMethod == "penalty")
    {
            // set initial quess for penalty factors
	    Eigen::MatrixXd tauInit(4,1);
            tauInit << 10, 10, 10, 10;
            //reduced.tauU = reduced.penalty_PPE_time(vel_now, tauInit);
            reduced.tauU = tauInit;
    }
   
    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (1); k++)
    {
	auto start_ROM = std::chrono::high_resolution_clock::now();
        //vel_now(0, 0) = par_on(0,k);
        reduced.solveOnline_PPE(vel_now, k);
        auto finish_ROM = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM = finish_ROM - start_ROM;
        std::cout << "elapsed_ROM: " << elapsed_ROM.count() << " seconds.";
        std::cout << std::endl;

        auto start_ROM_REC = std::chrono::high_resolution_clock::now();
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionPPE", 60);
        auto finish_ROM_REC = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_ROM_REC = finish_ROM_REC - start_ROM_REC;
        std::cout << "elapsed_ROM_REC: " << elapsed_ROM_REC.count() << " seconds.";
        std::cout << std::endl;
    }

/*   
    // Performing full order simulation for second parameter set 
    tutorialY HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.0001;
    HFonline2.mu_range(0, 1) = 0.0001;
    HFonline2.genEquiPar();
    HFonline2.inletIndex.resize(1, 2);
    HFonline2.inletIndex(0, 0) = 0;
    HFonline2.inletIndex(0, 1) = 0;
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 10;
    HFonline2.timeStep = 0.0005;
    HFonline2.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline2.onlineSolveFull(par_on, 1,
                              "./ITHACAoutput/HFonline2");



    // Performing full order simulation for second parameter set 
    tutorialY HFonline3(argc, argv);
    HFonline3.Pnumber = 1;
    HFonline3.Tnumber = 1;
    HFonline3.setParameters();
    HFonline3.mu_range(0, 0) = 0.0001;
    HFonline3.mu_range(0, 1) = 0.0001;
    HFonline3.genEquiPar();
    HFonline3.inletIndex.resize(1, 2);
    HFonline3.inletIndex(0, 0) = 0;
    HFonline3.inletIndex(0, 1) = 0;
    HFonline3.startTime = 0.0;
    HFonline3.finalTime = 10;
    HFonline3.timeStep = 0.0005;
    HFonline3.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline3.onlineSolveFull(par_on, 2,
                              "./ITHACAoutput/HFonline3"); */

 /*   // Reading high-fidelity solutions for the parameter set
    // for which the offline solve has been performed (skipping IC)
    example.onlineSolveRead("./ITHACAoutput/Offline/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline3/");

    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::error_listfields_min_IC(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixP = ITHACAutilities::error_listfields_min_IC(
                                         example.Pfield_on, reduced.PREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixP, "L2errorMatrixP", "eigen",
                               "./ITHACAoutput/l2error");




//Post-Process
    Eigen::MatrixXd PostP(example.Ufield_on.size(), 6); 

    for (label i = 0; i < example.Ufield_on.size(); i++)
    {
	PostP(i, 0) =  0.5*fvc::domainIntegrate(example.Ufield_on[i] & example.Ufield_on[i]).value();
	
	// min/max Outlet (patch 1) velocity
	PostP(i, 1) = max(example.Ufield_on[i].boundaryField()[1].component(0) );
        PostP(i, 2) = min(example.Ufield_on[i].boundaryField()[1].component(0) );

	PostP(i, 3) = 0.5*fvc::domainIntegrate(reduced.UREC[i] & reduced.UREC[i]).value();
	PostP(i, 4) = max(reduced.UREC[i].boundaryField()[1].component(0) );
        PostP(i, 5) = min(reduced.UREC[i].boundaryField()[1].component(0) );
	
    }

    ITHACAstream::exportMatrix(PostP, "PostP", "eigen", "./ITHACAoutput/PostProcess");

*/
exit(0); 
} 



  



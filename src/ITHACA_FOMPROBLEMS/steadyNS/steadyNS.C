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

\*---------------------------------------------------------------------------*/


/// \file
/// Source file of the steadyNS class.

#include "steadyNS.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
steadyNS::steadyNS() {}
steadyNS::steadyNS(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();

 ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
#include "createFields.H"
#include "createFvOptions.H"
    turbulence->validate();
   
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none", 
             "The BC method must be set to lift or penalty in ITHACAdict");
    ExplicitMethod = ITHACAdict->lookupOrDefault<word>("ExplicitMethod", "A");
    M_Assert(ExplicitMethod == "Ales" || ExplicitMethod == "A" || ExplicitMethod == "B", 
             "The Explicit method must be set to Ales or A or B in ITHACAdict");
    PoissonMethod = ITHACAdict->lookupOrDefault<word>("PoissonMethod", "A");
    M_Assert(PoissonMethod == "FOM" || PoissonMethod == "ROM", 
             "The Explicit method must be set to FOM or ROM in ITHACAdict");
    para = new ITHACAparameters;
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to perform a truthSolve
void steadyNS::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolve.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U);
    Pfield.append(p);
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (int i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}

// Method to solve the supremizer problem
void steadyNS::solvesupremizer(word type)
{
    M_Assert(type == "modes"
             || type == "snapshots",
             "You must specify the variable type with either snapshots or modes");
    PtrList<volScalarField> P_sup;

    if (type == "modes")
    {
        P_sup = Pmodes;
    }
    else if (type == "snapshots")
    {
        P_sup = Pfield;
    }

    if (supex == 1)
    {
        volVectorField U = _U();
        volVectorField Usup
        (
            IOobject
            (
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
        );

        if (type == "snapshots")
        {
            ITHACAstream::read_fields(supfield, Usup, "./ITHACAoutput/supfield/");
        }
        else
        {
            ITHACAstream::read_fields(supmodes, Usup, "./ITHACAoutput/supremizer/");
        }
    }
    else
    {
        volVectorField U = _U();
        volVectorField Usup
        (
            IOobject
            (
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
        );
        dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1)
        );
        Vector<double> v(0, 0, 0);

        for (label i = 0; i < Usup.boundaryField().size(); i++)
        {
            if (Usup.boundaryField()[i].type() != "processor")
            {
                ITHACAutilities::changeBCtype(Usup, "fixedValue", i);
                assignBC(Usup, i, v);
                assignIF(Usup, v);
            }
        }

        if (type == "snapshots")
        {
            for (label i = 0; i < P_sup.size(); i++)
            {
                fvVectorMatrix u_sup_eqn
                (
                    - fvm::laplacian(nu_fake, Usup)
                );
                solve
                (
                    u_sup_eqn == fvc::grad(P_sup[i])
                );
                supfield.append(Usup);
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfield/");
            }

            ITHACAutilities::createSymLink("./ITHACAoutput/supfield");
        }
        else
        {
            for (label i = 0; i < Pmodes.size(); i++)
            {
                fvVectorMatrix u_sup_eqn
                (
                    - fvm::laplacian(nu_fake, Usup)
                );
                solve
                (
                    u_sup_eqn == fvc::grad(Pmodes[i])
                );
                supmodes.append(Usup);
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supremizer/");
            }

            ITHACAutilities::createSymLink("./ITHACAoutput/supremizer");
        }
    }
}

// Method to compute the lifting function
void steadyNS::liftSolve()
{

    for (label k = 0; k < inletIndex.rows(); k++)
    {
Info << "bug1a" << endl;
        Time& runTime = _runTime();
Info << "bug1b" << endl;
	volScalarField p = _p();
Info << "bug1c" << endl;
        volVectorField U = _U();
        surfaceScalarField& phi = _phi();
Info << "bug2" << endl;
        fvMesh& mesh = _mesh();
Info << "bug3" << endl;
        IOMRFZoneList& MRF = _MRF();
Info << "bug4" << endl;
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
Info << "bug5" << endl;
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
Info << "bug6" << endl;
        pisoControl potentialFlow(mesh, "potentialFlow");
Info << "bug7" << endl;
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
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
Info << "bugA" << endl;
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
Info << "bugB" << endl;
       /* setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );*/
Info << "bugC" << endl;
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);
Info << "bugD" << endl;
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
Info << "bug0" << endl;
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

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void steadyNS::projectPPE(fileName folder, label NU, label NP, label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        }
        else
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word D_str = "D_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + D_str))
        {
            ITHACAstream::ReadDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/", D_str);
        }
        else
        {
            D_matrix = laplacian_pressure(NPmodes);
        }

        word bc3_str = "BC3_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc3_str))
        {
            ITHACAstream::ReadDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/", bc3_str);
        }
        else
        {
            BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        }

        word bc4_str = "BC4_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc4_str))
        {
            ITHACAstream::ReadDenseMatrix(BC4_matrix, "./ITHACAoutput/Matrices/", bc4_str);
        }
        else
        {
            BC4_matrix = pressure_BC4(NUmodes, NPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        word G_str = "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + G_str))
        {
            ITHACAstream::ReadDenseTensor(gTensor, "./ITHACAoutput/Matrices/", G_str);
        }
        else
        {
            gTensor = divMomentum(NUmodes, NPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        D_matrix = laplacian_pressure(NPmodes);
        gTensor = divMomentum(NUmodes, NPmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        BC4_matrix = pressure_BC4(NUmodes, NPmodes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "eigen",
                                   "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
    }
}

void steadyNS::projectSUP(fileName folder, label NU, label NP, label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
   


    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(supmodes[k]);
        }
    }
//Change back
    if (bcMethod == "none")
    {
	for (label j = 0; j < NUmodes; j++)
        {

	    	forAll( L_U_SUPmodes[j].boundaryFieldRef()[0], l)
	    	{
			L_U_SUPmodes[j].boundaryFieldRef()[0][l].component(vector::X) = 0;
			L_U_SUPmodes[j].boundaryFieldRef()[0][l].component(vector::Y) = 0;
	        	L_U_SUPmodes[j].boundaryFieldRef()[0][l].component(vector::Z) = 0;
	    	}
	}	    
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word P_str = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        }

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        }
        else
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }

	C_matrix  = convective_term(NUmodes, NPmodes, NSUPmodes);
	BC_matrix = boundary_term(NUmodes, NSUPmodes);
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
	//C_matrix  = convective_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
	BC_matrix = boundary_term(NUmodes, NSUPmodes);
	if (ExplicitMethod == "Ales")
        {
	     Bbc_matrix = diffusive_term_BC(NUmodes, NPmodes, NSUPmodes);
	    RedLinSysP = pressure_gradient_term_linsys(NPmodes);
	    RedLinSysPDiff = pressure_gradient_term_linsys_diff(NPmodes);
	    RedLinSysPConv = pressure_gradient_term_linsys_conv(NPmodes);
	    DF_matrix = laplacian_pressure_FOM(NUmodes, NPmodes);
	    Cf_tensor = convective_flux_tens(NUmodes, NPmodes, NSUPmodes);
	    //BC_matrix_PPE = boundary_term_PPE(NPmodes);
	    //Pf_matrix = divergence_flux_term(NUmodes, NPmodes, NSUPmodes);
	}
	else if (ExplicitMethod == "A")
        {
	    BC_matrix_PPE = boundary_term_PPE(NPmodes);
	   // BC_matrix_PPE2 = boundary_term_PPE2(NPmodes);
	    Cf_tensor = convective_flux_tens(NUmodes, NPmodes, NSUPmodes);
	    //Pf_matrix = divergence_flux_term(NUmodes, NPmodes, NSUPmodes);
	    DF_matrix = laplacian_pressure_FOM(NUmodes, NPmodes);
	    RedLinSysP = pressure_gradient_term_linsys(NPmodes);
	    RedLinSysPDiff = pressure_gradient_term_linsys_diff(NPmodes);
	    RedLinSysPConv = pressure_gradient_term_linsys_conv(NPmodes);
	    Kf_matrix = pressure_flux_term(NUmodes, NPmodes, NSUPmodes);
	    I_matrix = interpolation_term(NUmodes,NPmodes, NSUPmodes);
	    Mf_matrix = mass_matrix_flux(NUmodes,NPmodes, NSUPmodes);
	    ID_matrix = interpolation_term_diffusion(NUmodes,NPmodes, NSUPmodes);
	    Ci_matrix = convective_flux_int_mat(NUmodes, NPmodes, NSUPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python",
                                   "./ITHACAoutput/Matrices/C");
    }
}

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::diffusive_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label Bsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd B_matrix;
    B_matrix.resize(Bsize, Bsize);

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {

            B_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), L_U_SUPmodes[j]))).value();
        }
    }



   /* PtrList<volVectorField> diffU;
    for (int i = 0; i < Bsize; i++)
    {
         diffU.append(fvc::laplacian(dimensionedScalar("1", dimless, 1), L_U_SUPmodes[i]));
    }

    B_matrix = Modes<vector>::project(L_U_SUPmodes,diffU);*/

    if (Pstream::parRun())
    {
        reduce(B_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/",
                                  "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(B_matrix, "Br_", "eigen",
                               "./ITHACAoutput/Matrices/");

    return B_matrix;
}


// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::diffusive_term_BC(label NUmodes, label NPmodes,
        label NSUPmodes)
{

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	volVectorField U0 = _U();
    	Vector<double> inl(0, 0, 0);
    	assignIF(U0, inl);
 	

    	Eigen::MatrixXd Bbc_matrix(NPmodes,NPmodes);

   	// Project everything
    	for (label i = 0; i < NPmodes; i++)
    	{
	    for (label j = 0; j < NPmodes; j++)
    	    {

		volVectorField L_U_SUPmodesj(_Ub);
		L_U_SUPmodesj = L_U_SUPmodes[j];

		volVectorField UBCU(_U0);
	    	UBCU = dt_fake*(-fvc::div(fvc::flux(U0),L_U_SUPmodesj)-fvc::div(fvc::flux(L_U_SUPmodesj),U0));
		//volScalarField UBCP = (1/dt_fake)*fvc::div(UBCU);
		Bbc_matrix(i,j) = fvc::domainIntegrate(Pmodes[i] *  fvc::div(UBCU)).value();
	    }

    	}


    ITHACAstream::exportMatrix(Bbc_matrix, "Bbcr_", "eigen",
                               "./ITHACAoutput/Matrices/");


    return Bbc_matrix;
}

Eigen::MatrixXd steadyNS::mass_matrix_flux(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    Eigen::MatrixXd Mf_matrix;
    Mf_matrix.resize(NUmodes,NUmodes);

    for (label i = 0; i < NUmodes; i++)
    {
	volVectorField A = fvc::reconstruct(Phimodes[i]);
	for (label j = 0; j < NUmodes; j++)
        {
	    volVectorField B = fvc::reconstruct(Phimodes[j]);
	    Mf_matrix(i, j) = fvc::domainIntegrate(A & B).value();

	}
    }

   /* PtrList<volVectorField> A;
    for (int i = 0; i < Isize; i++)
    {
         A.append(fvc::reconstruct(Phimodes[i]));
    }

    Mf_matrix = A.project(A);*/


    ITHACAstream::SaveDenseMatrix(Mf_matrix, "./ITHACAoutput/Matrices/",
                                  "Mf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(Mf_matrix, "Mphir_", "eigen",
                               "./ITHACAoutput/Matrices/");

    return Mf_matrix;
}

Eigen::MatrixXd steadyNS::interpolation_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{

    label Isize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd I_matrix;
    I_matrix.resize(NUmodes, Isize);

    volVectorField U0 = _U();
    Vector<double> inl(0, 0, 0);
    assignIF(U0, inl);
    volVectorField CoeffA("CoeffA", U0);
    volVectorField CoeffB("CoeffB", U0);
    for (label l = 0; l < NUmodes; l++)
    {
	CoeffA = fvc::reconstruct(Phimodes[l]);
	for (label j = 0; j < Isize; j++)
        {
	    CoeffB = fvc::reconstruct(fvc::flux(L_U_SUPmodes[j]));
	    I_matrix(l,j) = fvc::domainIntegrate(CoeffA &  CoeffB).value();
	}
    }

    /*PtrList<volVectorField> A;
    PtrList<volVectorField> B;
    for (int i = 0; i < Isize; i++)
    {
         B.append(fvc::reconstruct(Phimodes[i]));
	 A.append(fvc::reconstruct(fvc::flux(Umodes[i])));
    }

    I_matrix = A.project(B);*/

    if (Pstream::parRun())
    {
        reduce(I_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(I_matrix, "./ITHACAoutput/Matrices/",
                                  "I_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(I_matrix, "Ir_", "eigen",
                               "./ITHACAoutput/Matrices/");

    return I_matrix;
}

Eigen::MatrixXd steadyNS::interpolation_term_diffusion(label NUmodes, label NPmodes,
        label NSUPmodes)
{

    dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );

    label Isize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd ID_matrix;
    ID_matrix.resize(Isize, NUmodes);
    volVectorField U0 = _U();
    Vector<double> inl(0, 0, 0);
    assignIF(U0, inl);

    volVectorField CoeffC("CoeffC", _Ub);
    surfaceScalarField phi_4("Phi_4", _phi);
    for (label l = 0; l < NUmodes; l++)
    {
	for (label j = 0; j < Isize; j++)
        {
	   CoeffC =  Umodes[j];
	   // CoeffC = L_U_SUPmodes[j];
	    phi_4 = dt_fake*nu_fake*fvc::flux(fvc::laplacian(dimensionedScalar("1", dimless, 1)
                   , CoeffC));
	    volVectorField CoeffB = fvc::reconstruct(phi_4);
	    volVectorField CoeffA = fvc::reconstruct(Phimodes[l]);
	    ID_matrix(l,j) = fvc::domainIntegrate(CoeffA &  CoeffB).value();
	}
    }

  /*  for (label l = 0; l < Isize; l++)
    {
	for (label j = 0; j < Isize; j++)
        {
	   //volVectorField CoeffA = fvc::reconstruct(fvc::flux(dt_fake*fvc::laplacian(
             //       nu_fake, Umodes[j])));
	    CoeffA = fvc::reconstruct(fvc::flux(dt_fake*fvc::laplacian(
                   nu_fake, Umodes[j])));
	   //  CoeffA = dt_fake*fvc::laplacian(
             //       nu_fake, L_U_SUPmodes[j]);
	    CoeffB = fvc::reconstruct(Phimodes[l]);
	    ID_matrix(l,j) = fvc::domainIntegrate(CoeffA &  CoeffB).value();
	}
    }*/

    /*PtrList<volVectorField> A;
    PtrList<volVectorField> B;
    for (int i = 0; i < Isize; i++)
    {
         B.append(fvc::reconstruct(Phimodes[i]));
	 A.append(fvc::reconstruct(fvc::flux(dt_fake*fvc::laplacian(
                    nu_fake, Umodes[i]))));
    }

    ID_matrix = A.project(B);*/

    if (Pstream::parRun())
    {
        reduce(ID_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(ID_matrix, "./ITHACAoutput/Matrices/",
                                  "ID_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(ID_matrix, "IDr_", "eigen",
                               "./ITHACAoutput/Matrices/");

    return ID_matrix;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys(label NPmodes)
{
	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );
  
	//volVectorField U(_U);

  
 	/*Vector<double> inl(0, 0, 0);
    	volVectorField U0 = _U();
    	assignIF(U0, inl);*/

	 volVectorField U0("U0",_U());
            Vector<double> inl(0, 0, 0);
           assignIF(U0, inl);
	    ITHACAutilities::changeBCtype( U0,"fixedValue",1);
    	    assignBC( U0,1,inl);


	//volVectorField U1aux(_U);
	//U1aux = U0 + dt_fake * (fvc::laplacian(nu_fake,U0)-fvc::div(fvc::flux(U0),U0));
	//volVectorField U1aux(_Ub);
	//U1aux = U0;
	volScalarField p(_p);
	
	fvScalarMatrix pEqn
        (
	    // fvm::laplacian(p) == (1/dt_fake)*fvc::div(U0)//  closed== (1/dt_fake)*fvc::div(U1aux)
	     fvm::laplacian(p) == (1/dt_fake)*fvc::div(U0)//  open== (1/dt_fake)*fvc::div(U1aux)
        );
        pEqn.setReference(0, 0);
        RedLinSysP = Pmodes.project(pEqn,NPmodes);

	return RedLinSysP;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys_conv(label NPmodes)
{
	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );
 
	volScalarField p(_p);
	/*volVectorField U0(_U);
 	Vector<double> inl(0, 0, 0);
 	assignIF(U0, inl);
	volVectorField Caux(_U0); //closed*/

	 volVectorField U0("U0",_U());
            Vector<double> inl(0, 0, 0);
           assignIF(U0, inl);
	    ITHACAutilities::changeBCtype( U0,"fixedValue",1);
    	    assignBC( U0,1,inl);
	volVectorField Caux(_Ub); 
	Caux = dt_fake * (-fvc::div(fvc::flux(U0),U0));
	
	fvScalarMatrix pEqn
        (
	     fvm::laplacian(p) == (1/dt_fake)*fvc::div(Caux)
        );
        pEqn.setReference(0, 0);
        RedLinSysPConv = Pmodes.project(pEqn,NPmodes);

	return RedLinSysPConv;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys_diff(label NPmodes)
{
        dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );
	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );
  
	/*//volVectorField U(_U);
	volVectorField U0(_U);
 	Vector<double> inl(0, 0, 0);
 	assignIF(U0, inl);
	volVectorField Daux(_U0);
	//volVectorField Daux(_Ub);//closed*/

	 volVectorField U0("U0",_U());
            Vector<double> inl(0, 0, 0);
           assignIF(U0, inl);
	    ITHACAutilities::changeBCtype( U0,"fixedValue",1);
    	    assignBC( U0,1,inl);
	volVectorField Daux(_Ub); 
	Daux = dt_fake * fvc::laplacian(nu_fake,U0);

	volScalarField p(_p);
	
	fvScalarMatrix pEqn
        (
             fvm::laplacian(p) == (1/dt_fake)*fvc::div(Daux)

        );
        pEqn.setReference(0, 0);
        RedLinSysPDiff = Pmodes.project(pEqn,NPmodes);

	return RedLinSysPDiff;
}

Eigen::MatrixXd steadyNS::pressure_gradient_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label K1size = NUmodes + NSUPmodes + liftfield.size();
    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
	
            K_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::grad(
                    Pmodes[j])).value();
        }
    }


  /*  PtrList<volVectorField> A;
    PtrList<volVectorField> B;
    for (int i = 0; i < K1size; i++)
    {
	 A.append(L_U_SUPmodes[i]);
    }

    for (int i = 0; i < K2size; i++)
    {
         B.append(fvc::grad(Pmodes[i]));
    }

    K_matrix = A.project(B);*/

    if (Pstream::parRun())
    {
        reduce(K_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/",
                                  "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));

    ITHACAstream::exportMatrix(K_matrix, "Kr_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return K_matrix;
}

Eigen::MatrixXd steadyNS::pressure_flux_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    //label Kf1size = NUmodes + NSUPmodes + liftfield.size();
    label Kf1size = NUmodes ;
    label Kf2size = NPmodes;
    Eigen::MatrixXd Kf_matrix(Kf1size, Kf2size);

     dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );

    for (label l = 0; l < Kf1size; l++)
    {
	for (label j = 0; j < Kf2size; j++)
        {
	    volVectorField CoeffA = fvc::reconstruct(dt_fake*fvc::snGrad(Pmodes[j])*mag(Pmodes[j].mesh().magSf()));
	    volVectorField CoeffB = fvc::reconstruct(Phimodes[l]);
	    Kf_matrix(l,j) = fvc::domainIntegrate(CoeffA &  CoeffB).value();
	}
    }

   /* PtrList<volVectorField> A;
    PtrList<volVectorField> B;
    for (int i = 0; i < Kf2size; i++)
    {
	 A.append(fvc::reconstruct(dt_fake*fvc::snGrad(Pmodes[i])*mag(Pmodes[i].mesh().magSf())));
    }

    for (int i = 0; i < Kf1size; i++)
    {
         B.append(fvc::reconstruct(Phimodes[i]));
    }

    Kf_matrix = A.project(B);*/


    if (Pstream::parRun())
    {
        reduce(Kf_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(Kf_matrix, "./ITHACAoutput/Matrices/",
                                  "Kf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));

    ITHACAstream::exportMatrix(Kf_matrix, "Kfr_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return Kf_matrix;
}

List <Eigen::MatrixXd> steadyNS::convective_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{

    Info << "C_matrix" << endl;
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> C_matrix;
    C_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        C_matrix[j].resize(Csize, Csize);
    }

    for (label i = 0; i < Csize; i++)
    {
        for (label j = 0; j < Csize; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                C_matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf() ,
                                        L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        for (int i = 0; i < Csize; i++)
        {
            List<double> vec(C_matrix[i].data(), C_matrix[i].data() + C_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(C_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(C_matrix, "Cr", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(C_matrix, "Cr", "eigen", "./ITHACAoutput/Matrices/Cr");
    return C_matrix;
}

Eigen::Tensor<double, 3> steadyNS::convective_term_tens(label NUmodes,
        label NPmodes,
        label NSUPmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> C_tensor;
    

	    if (ExplicitMethod == "Ales")
	    {
		C_tensor.resize(Csize, Csize, Csize);

		 for (label i = 0; i < Csize; i++)
    		{
       
			for (label j = 0; j < Csize; j++)
        		{
                		for (label k = 0; k < Csize; k++)
                		{
		
                     			C_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_U_SUPmodes[k])).value();
				}
			}
		}
	    }
	    else if (ExplicitMethod == "A")
	    {

		C_tensor.resize(Csize, NUmodes, Csize);
	 	for (label i = 0; i < Csize; i++)
    		{
       
			for (label j = 0; j < NUmodes; j++)
        		{
				for (label k = 0; k < Csize; k++)
                		{
	             		C_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                        Phimodes[j],
                                        L_U_SUPmodes[k])).value();
				}
	    		}
	
           	}
           }


    if (Pstream::parRun())
    {
        reduce(C_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(C_tensor, "./ITHACAoutput/Matrices/",
                                  "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return C_tensor;
}

Eigen::Tensor<double, 3> steadyNS::convective_flux_tens(label NUmodes,
        label NPmodes, label NSUPmodes)
{

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );


    label Csize1 = NUmodes + NSUPmodes + liftfield.size();
    label Csize2 = NPmodes;
    Eigen::Tensor<double, 3> Cf_tensor;


if (ExplicitMethod == "Ales")
{
    Cf_tensor.resize(Csize2, Csize1, Csize1);
    for (label i = 0; i < Csize2; i++)
    {
        for (label j = 0; j < Csize1; j++)
        {

            for (label k = 0; k < Csize1; k++)
            {
	     
	        volVectorField L_U_SUPmodesaux(_Ub);
		//volVectorField L_U_SUPmodesaux(_U);
		
                    L_U_SUPmodesaux = dt_fake*(fvc::div(
                                        fvc::flux(L_U_SUPmodes[j] ) ,
                                       L_U_SUPmodes[k]));

		Cf_tensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * fvc::div(L_U_SUPmodesaux)).value();
	    }
	}
    }
}
else if (ExplicitMethod == "A")
{
    Cf_tensor.resize(Csize2, NUmodes, Csize1);
    for (label i = 0; i < Csize2; i++)
    {
        for (label j = 0; j < NUmodes; j++)
        {

            for (label k = 0; k < Csize1; k++)
            {
	     
	        volVectorField L_U_SUPmodesaux(_Ub);
		
                   L_U_SUPmodesaux = dt_fake*fvc::div(
                                        Phimodes[j],
                                        L_U_SUPmodes[k]);
		Cf_tensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * fvc::div(L_U_SUPmodesaux)).value();
	    }
	}
    }
}


    if (Pstream::parRun())
    {
        reduce(Cf_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(Cf_tensor, "./ITHACAoutput/Matrices/",
                                  "Cf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return Cf_tensor;
}
List <Eigen::MatrixXd> steadyNS::convective_flux_int_mat(label NUmodes, label NPmodes,
        label NSUPmodes)
{

    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> Ci_matrix;
    Ci_matrix.setSize(NUmodes);

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );


    for (label j = 0; j < NUmodes; j++)
    {
        Ci_matrix[j].resize(NUmodes, Csize);
    }

    volVectorField U0 = _U();
    Vector<double> inl(0, 0, 0);
    assignIF(U0, inl);

    volVectorField CoeffD("CoeffD", _U);
    surfaceScalarField phi_4("Phi_4", _phi);
    surfaceScalarField CoeffC("CoeffC", _phi);
    for (label i = 0; i < NUmodes; i++)
    {
        for (label j = 0; j < NUmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
		CoeffC = Phimodes[i];
		CoeffD = L_U_SUPmodes[k];
		//CoeffD = Umodes[k];
		phi_4 = dt_fake*fvc::flux(fvc::div(CoeffC,CoeffD));
	        volVectorField CoeffA = fvc::reconstruct(phi_4);
	      //volVectorField CoeffA = dt_fake*fvc::div(Phimodes[i],L_U_SUPmodes[k]);
	        volVectorField CoeffB = fvc::reconstruct(Phimodes[j]);
	        Ci_matrix[i](j, k) = fvc::domainIntegrate(CoeffB &  CoeffA).value();
            }
        }
    }


    // Export the matrix
    ITHACAstream::exportMatrix(Ci_matrix, "Cir", "eigen", "./ITHACAoutput/Matrices/Cir");
    return Ci_matrix;
}

/*Eigen::Tensor<double, 3> steadyNS::convective_flux_int_tens(label NUmodes,
        label NPmodes, label NSUPmodes)
{
    label Csize1 = NUmodes + NSUPmodes + liftfield.size();
    label Csize2 = NUmodes;
    Eigen::Tensor<double, 3> Ci_tensor;
    Ci_tensor.resize(Csize2, Csize1, Csize1);

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );



for (label m = 0; m < Csize1; m++)
    {
	for (label l = 0; l < Csize1; l++)
        {
	    for (label j = 0; j < Csize1; j++)
            {
	        volVectorField CoeffA = fvc::reconstruct(fvc::flux(
			dt_fake*fvc::div(Phimodes[m],Umodes[j])));
	        volVectorField CoeffB = fvc::reconstruct(Phimodes[l]);
	        Ci_tensor(m, l, j) = fvc::domainIntegrate(CoeffA &  CoeffB).value();
	    }
	}
    }



    if (Pstream::parRun())
    {
        reduce(Ci_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(Ci_tensor, "./ITHACAoutput/Matrices/",
                                  "Ci_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return Ci_tensor;
}*/

Eigen::MatrixXd steadyNS::mass_term(label NUmodes, label NPmodes,
                                    label NSUPmodes)
{
    label Msize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd M_matrix(Msize, Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                                                  L_U_SUPmodes[j]).value();
        }
    }

   /* PtrList<volVectorField> A;

    for (int i = 0; i < Msize; i++)
    {
	 A.append(L_U_SUPmodes[i]);
    }

    M_matrix = A.project(A);*/

    if (Pstream::parRun())
    {
        reduce(M_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/",
                                  "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(M_matrix, "Mr_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return M_matrix;
}

Eigen::MatrixXd steadyNS::boundary_term(label NUmodes,
                                    label NSUPmodes)
{
    label BCsize = NUmodes + NSUPmodes + liftfield.size();
    word bw_str = "bw";
    word matname = "./ITHACAoutput/BCvector/" + bw_str + "0" + "_mat.txt";
    Eigen::MatrixXd bw = ITHACAstream::readMatrix(matname);

    Eigen::MatrixXd BC_matrix(BCsize,1);
    Eigen::VectorXd ModeVector;
 
    // Project everything
    for (label i = 0; i < BCsize; i++)
    {
	ModeVector = Foam2Eigen::field2Eigen(L_U_SUPmodes[i]);
	//Eigen::VectorXd ModeVector2 = Foam2Eigen::field2Eigen(L_U_SUPmodes[i].mesh().V());
	//Eigen::VectorXd ModeVector3 = ModeVector * ModeVector2(0);
	//BC_matrix(i,0) = (ModeVector3.dot(bw.col(0)));
	BC_matrix(i,0) = ModeVector.dot(bw.col(0));

    }

    ITHACAstream::SaveDenseMatrix(BC_matrix, "./ITHACAoutput/Matrices/",
                                  "BC_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(BC_matrix, "BCr_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return BC_matrix;
}

Eigen::MatrixXd steadyNS::boundary_term_PPE(label NPmodes)
{

dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );
    label BCsize = NPmodes;

	volVectorField U0 = _U();
    Vector<double> inl(0, 0, 0);
    assignIF(U0, inl);
   surfaceScalarField phi_3("Phi_3", _phi);
    Eigen::MatrixXd BC_matrix_PPE(BCsize,1);

    // Project everything
    for (label i = 0; i < BCsize; i++)
    {
	volVectorField CoeffB = fvc::reconstruct(Phimodes[i]);
	phi_3 = dt_fake*nu_fake*fvc::flux(fvc::laplacian(dimensionedScalar("1", dimless, 1),U0));
	volVectorField CoeffA = fvc::reconstruct(phi_3);
	
	BC_matrix_PPE(i,0) = fvc::domainIntegrate(CoeffA &  CoeffB).value();

    }


    ITHACAstream::exportMatrix(BC_matrix_PPE, "BCr_PPE_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return BC_matrix_PPE;
}


Eigen::MatrixXd steadyNS::boundary_term_PPE2(label NPmodes)
{
    label BCsize = NPmodes;
dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );
    Eigen::MatrixXd BC_matrix_PPE2(BCsize,1);
   /* Eigen::VectorXd ModeVector;

 
	volVectorField PhiRec2 = fvc::reconstruct(Phifield[0]);
    // Project everything
    for (label i = 0; i < BCsize; i++)
    {
	 volVectorField PhiRec3 = fvc::reconstruct(Phimodes[i]);
	 
         BC_matrix_PPE2(i,0) = fvc::domainIntegrate(PhiRec3 &
                                                  PhiRec2).value();
    }*/

	volVectorField U0 = _U();
    Vector<double> inl(0, 0, 0);
    assignIF(U0, inl);
   surfaceScalarField phi_3("Phi_3", _phi);
    Eigen::MatrixXd BC_matrix_PPE(BCsize,1);

    // Project everything
    for (label i = 0; i < BCsize; i++)
    {
	volVectorField CoeffB = fvc::reconstruct(Phimodes[i]);
	phi_3 = dt_fake*fvc::flux(fvc::div(fvc::flux(U0),U0));
	volVectorField CoeffA = fvc::reconstruct(phi_3);
	
	BC_matrix_PPE2(i,0) = fvc::domainIntegrate(CoeffA &  CoeffB).value();

    }

   /* Eigen::MatrixXd M_matrixphi;
    M_matrixphi.resize(BCsize,BCsize);

    Eigen::MatrixXd a;
    a.resize(BCsize,1);

    // Project everything
    for (label i = 0; i < BCsize; i++)
    {
	volVectorField A = fvc::reconstruct(Phimodes[i]);
        for (label j = 0; j < BCsize; j++)
        {
	    volVectorField B = fvc::reconstruct(Phimodes[j]);
	    M_matrixphi(i, j) = fvc::domainIntegrate(A & B).value();
        }

    }

    a = M_matrixphi.colPivHouseholderQr().solve(BC_matrix_PPE2.col(0));

     BC_matrix_PPE2.col(0) = a;

    ITHACAstream::SaveDenseMatrix(BC_matrix_PPE2, "./ITHACAoutput/Matrices/",
                                  "BC_PPE2_" );*/

    ITHACAstream::exportMatrix(BC_matrix_PPE2, "BCr_PPE2_", "eigen",
                               "./ITHACAoutput/Matrices/");
    return BC_matrix_PPE2;
}


// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::divergence_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label P1size = NPmodes;
    label P2size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd P_matrix(P1size, P2size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
	    volVectorField U2aux(_Ub);
	    U2aux = L_U_SUPmodes[j];
	    //   Vector<double> inl(0, 0, 0);
	   // ITHACAutilities::changeBCtype( U2aux,"fixedValue",1);
    	    //ITHACAutilities::assignBC( U2aux,1,inl);
	
	   
            P_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                  U2aux)).value();
        }
    }

  /*  PtrList<volScalarField> A;
    PtrList<volScalarField> B;
    for (int i = 0; i < P1size; i++)
    {
	 A.append(Pmodes[i]);
    }
   for (int i = 0; i < P2size; i++)
    {
	 B.append(fvc::div(L_U_SUPmodes[i]));
    }

    P_matrix = A.project(B);*/

    if (Pstream::parRun())
    {
        reduce(P_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/",
                                  "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));

    ITHACAstream::exportMatrix(P_matrix, "Pr_", "eigen", "./ITHACAoutput/Matrices/");
    return P_matrix;
}

Eigen::MatrixXd steadyNS::divergence_term_pressure(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label P1size = NPmodes;
    Eigen::MatrixXd Pp_matrix(P1size, P1size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {

	
        for (label j = 0; j < P1size; j++)
        {
	
            Pp_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                   fvc::grad(Pmodes[j]))).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(Pp_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(Pp_matrix, "./ITHACAoutput/Matrices/",
                                  "Pp_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));

    ITHACAstream::exportMatrix(Pp_matrix, "Ppr_", "eigen", "./ITHACAoutput/Matrices/");
    return Pp_matrix;
}

Eigen::MatrixXd steadyNS::divergence_flux_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label P1size = NPmodes;
    label P2size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd Pf_matrix(P1size, P2size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
	            Pf_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                    fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf())).value();
        }
    }


    if (Pstream::parRun())
    {
        reduce(Pf_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(Pf_matrix, "./ITHACAoutput/Matrices/",
                                  "Pf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));

    ITHACAstream::exportMatrix(Pf_matrix, "Pfr_", "eigen", "./ITHACAoutput/Matrices/");
    return Pf_matrix;
}

List <Eigen::MatrixXd> steadyNS::div_momentum(label NUmodes, label NPmodes)
{
    label G1size = NPmodes;
    label G2size = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> G_matrix;
    G_matrix.setSize(G1size);

    for (label j = 0; j < G1size; j++)
    {
        G_matrix[j].resize(G2size, G2size);
    }

    for (label i = 0; i < G1size; i++)
    {
        for (label j = 0; j < G2size; j++)
        {
            for (label k = 0; k < G2size; k++)
            {
                G_matrix[i](j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_U_SUPmodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        for (int i = 0; i < G2size; i++)
        {
            List<double> vec(G_matrix[i].data(), G_matrix[i].data() + G_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(G_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(G_matrix, "G", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(G_matrix, "G", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(G_matrix, "G", "eigen", "./ITHACAoutput/Matrices/G");
    return G_matrix;
}

Eigen::Tensor<double, 3> steadyNS::divMomentum(label NUmodes, label NPmodes)
{
    label g1Size = NPmodes;
    label g2Size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> gTensor;
    gTensor.resize(g1Size, g2Size, g2Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            for (label k = 0; k < g2Size; k++)
            {
                gTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                        L_U_SUPmodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(gTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(gTensor, "./ITHACAoutput/Matrices/",
                                  "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return gTensor;
}

Eigen::Tensor<double, 3> steadyNS::divMomentumConvection(label NUmodes, label NPmodes)
{
    label g1Size = NPmodes;
    label g2Size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> gTensor;
    gTensor.resize(g1Size, g2Size, g2Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            for (label k = 0; k < g2Size; k++)
            {
                gFTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * 
				    fvc::div(
				    fvc::div(
                        (fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf()),
                        L_U_SUPmodes[k])
				    )
				    ).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(gFTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(gFTensor, "./ITHACAoutput/Matrices/",
                                  "GF_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return gFTensor;
}

Eigen::MatrixXd steadyNS::laplacian_pressure(label NPmodes)
{
    label Dsize = NPmodes;
    Eigen::MatrixXd D_matrix(Dsize, Dsize);

    // Project everything
    for (label i = 0; i < Dsize; i++)
    {
        for (label j = 0; j < Dsize; j++)
        {
            D_matrix(i, j) = fvc::domainIntegrate(fvc::grad(Pmodes[i])&fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(D_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/",
                                  "D_" + name(NPmodes));

     ITHACAstream::exportMatrix(D_matrix, "Dr_", "eigen", "./ITHACAoutput/Matrices/");

    return D_matrix;
}




Eigen::MatrixXd steadyNS::laplacian_pressure_FOM(label NUmodes, label NPmodes)
{
    label Dsize1 = NPmodes;
    label Dsize2 = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd DF_matrix(Dsize1, Dsize2);

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );

	dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1.0)
        );

  

    // Project everything
    for (label i = 0; i < Dsize1; i++)
    {
        for (label j = 0; j < Dsize2; j++)
        {

	    volVectorField L_U_SUPmodesaux(_Ub); //closed flow
	    //volVectorField L_U_SUPmodesaux(_U);
	    

	    L_U_SUPmodesaux = dt_fake*fvc::laplacian(
                    nu_fake, L_U_SUPmodes[j]);


            DF_matrix(i, j) = fvc::domainIntegrate(Pmodes[i]*fvc::div(L_U_SUPmodesaux)).value();
        }
    }

   /* PtrList<volScalarField> A;
    PtrList<volScalarField> B;
    for (int i = 0; i < Dsize1; i++)
    {
	 A.append(Pmodes[i]);
    }
   for (int i = 0; i < Dsize2; i++)
    {

	 B.append(dt_fake*fvc::laplacian(nu_fake, L_U_SUPmodes[i]));
    }

    DF_matrix = A.project(B);*/

    if (Pstream::parRun())
    {
        reduce(DF_matrix, sumOp<Eigen::MatrixXd>());
    }

  ITHACAstream::exportMatrix(DF_matrix, "DFr_", "eigen", "./ITHACAoutput/Matrices/");

    ITHACAstream::SaveDenseMatrix(DF_matrix, "./ITHACAoutput/Matrices/",
                                  "DF_" + name(NPmodes));
    return DF_matrix;
}



/*Eigen::MatrixXd steadyNS::laplacian_pressure_bc(label NUmodes, label NPmodes,
        label NSUPmodes)
{

	Eigen::MatrixXd A;
    	Eigen::VectorXd b;

	volVectorField Diffbc0 = fvc::laplacian(dimensionedScalar("1", dimless, 1), L_U_SUPmodes[0]);
	volScalarField DiffbcScalar = fvc::div(Diffbc0);
	fvScalarMatrix pEqn0
    	(
		fvc::div(Diffbc0)
    	);
	Foam2Eigen::fvMatrix2Eigen(pEqn0, A, b);

	Eigen::MatrixXd Mode_matrix; Mode_matrix.resize(b.size(), NPmodes);


	

 for (label i = 0; i < NUmodes; i++)
    {
	A = A * 0;
	b= b * 0;

	volVectorField Diffbc1 = fvc::laplacian(dimensionedScalar("1", dimless, 1), L_U_SUPmodes[i]);
	fvScalarMatrix PEqn
    	(
		fvc::div(Diffbc1)
    	);
	Foam2Eigen::fvMatrix2Eigen(PEqn, A, b);
	Mode_matrix.col(i) = b;
    
    }
	Eigen::MatrixXd BCrhs_matrix = Foam2Eigen::PtrList2Eigen(Pmodes);
	DFbc_matrix = Mode_matrix.transpose()*(BCrhs_matrix);


    ITHACAstream::exportMatrix(DFbc_matrix, "DFbcr_", "eigen",
                               "./ITHACAoutput/Matrices/");

    return DFbc_matrix;
}*/


Eigen::MatrixXd steadyNS::pressure_BC1(label NUmodes, label NPmodes)
{
    label P_BC1size = NPmodes;
    label P_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
    fvMesh& mesh = _mesh();

    for (label i = 0; i < P_BC1size; i++)
    {
        for (label j = 0; j < P_BC2size; j++)
        {
            surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(
                                        L_U_SUPmodes[j]))&mesh.Sf())*fvc::interpolate(Pmodes[i]));
            double s = 0;

            for (label k = 0; k < lpl.boundaryField().size(); k++)
            {
                s += gSum(lpl.boundaryField()[k]);
            }

            BC1_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC1_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/",
                                  "BC1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC1_matrix;
}


List <Eigen::MatrixXd> steadyNS::pressure_BC2(label NUmodes, label NPmodes)
{
    label P2_BC1size = NPmodes;
    label P2_BC2size = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> BC2_matrix;
    fvMesh& mesh = _mesh();
    BC2_matrix.setSize(P2_BC1size);

    for (label j = 0; j < P2_BC1size; j++)
    {
        BC2_matrix[j].resize(P2_BC2size, P2_BC2size);
    }

    for (label i = 0; i < P2_BC1size; i++)
    {
        for (label j = 0; j < P2_BC2size; j++)
        {
            for (label k = 0; k < P2_BC2size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                             L_U_SUPmodes[j]) & mesh.Sf(),
                                         L_U_SUPmodes[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    s += gSum(div_m.boundaryField()[k]);
                }

                BC2_matrix[i](j, k) = s;
            }
        }
    }

    if (Pstream::parRun())
    {
        for (int i = 0; i < P2_BC1size; i++)
        {
            List<double> vec(BC2_matrix[i].data(),
                             BC2_matrix[i].data() + BC2_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(BC2_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    return BC2_matrix;
}

Eigen::Tensor<double, 3 > steadyNS::pressureBC2(label NUmodes, label NPmodes)
{
    label pressureBC1Size = NPmodes;
    label pressureBC2Size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3 > bc2Tensor;
    fvMesh& mesh = _mesh();
    bc2Tensor.resize(pressureBC1Size, pressureBC2Size, pressureBC2Size);

    for (label i = 0; i < pressureBC1Size; i++)
    {
        for (label j = 0; j < pressureBC2Size; j++)
        {
            for (label k = 0; k < pressureBC2Size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                             L_U_SUPmodes[j]) & mesh.Sf(),
                                         L_U_SUPmodes[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    s += gSum(div_m.boundaryField()[k]);
                }

                bc2Tensor(i, j, k) = s;
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(bc2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/",
                                  "BC2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return bc2Tensor;
}

Eigen::MatrixXd steadyNS::pressure_BC3(label NUmodes, label NPmodes)
{
    label P3_BC1size = NPmodes;
    label P3_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
    //fvMesh& mesh = _mesh();
    surfaceVectorField n(Pmodes[0].mesh().Sf() / Pmodes[0].mesh().magSf());

    for (label i = 0; i < P3_BC1size; i++)
    {
        for (label j = 0; j < P3_BC2size; j++)
        {
            surfaceVectorField BC3 = fvc::interpolate(fvc::curl(L_U_SUPmodes[j]));
            surfaceVectorField BC4 = n ^ fvc::interpolate(fvc::grad(Pmodes[i]));
            surfaceScalarField BC5 = (BC3 & BC4) * BC3.mesh().magSf();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                s += gSum(BC5.boundaryField()[k]);
            }

            BC3_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC3_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/",
                                  "BC3_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC3_matrix;
}

Eigen::MatrixXd steadyNS::pressure_BC4(label NUmodes, label NPmodes)
{
    label P4_BC1size = NPmodes;
    label P4_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC4_matrix(P4_BC1size, P4_BC2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P4_BC1size; i++)
    {
        for (label j = 0; j < P4_BC2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Pmodes[i]);
            surfaceScalarField BC4 = n & fvc::interpolate(L_U_SUPmodes[j]);
            surfaceScalarField BC5 = (BC3 * BC4) * mesh.magSf();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                s += gSum(BC5.boundaryField()[k]);
            }

            BC4_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC4_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC4_matrix, "./ITHACAoutput/Matrices/",
                                  "BC4_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC4_matrix;
}

List< Eigen::MatrixXd > steadyNS::bcVelocityVec(label NUmodes,
        label NSUPmodes)
{
    label BCsize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > bcVelVec(inletIndex.rows());

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelVec[j].resize(BCsize, 1);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < BCsize; i++)
        {
            bcVelVec[k](i, 0) = gSum(L_U_SUPmodes[i].boundaryField()[BCind].component(
                                         BCcomp));
        }
    }

    ITHACAstream::exportMatrix(bcVelVec, "bcVelVec", "eigen",
                               "./ITHACAoutput/Matrices/bcVelVec");
    return bcVelVec;
}

List< Eigen::MatrixXd > steadyNS::bcVelocityMat(label NUmodes,
        label NSUPmodes)
{
    label BCsize = NUmodes + NSUPmodes + liftfield.size() ;
    label BCUsize = inletIndex.rows();
    List < Eigen::MatrixXd > bcVelMat(BCUsize);

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelMat[j].resize(BCsize, BCsize);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < BCsize; i++)
        {
            for (label j = 0; j < BCsize; j++)
            {
                bcVelMat[k](i, j) = gSum(L_U_SUPmodes[i].boundaryField()[BCind].component(
                                             BCcomp) *
                                         L_U_SUPmodes[j].boundaryField()[BCind].component(BCcomp));
            }
        }
    }

    ITHACAstream::exportMatrix(bcVelMat, "bcVelMat", "eigen",
                               "./ITHACAoutput/Matrices/bcVelMat");
    return bcVelMat;
}

void steadyNS::change_viscosity(double mu)
{
    const volScalarField& nu =  _laminarTransport().nu();
    volScalarField& ciao = const_cast<volScalarField&>(nu);
    this->assignIF(ciao, mu);

    for (int i = 0; i < ciao.boundaryFieldRef().size(); i++)
    {
        this->assignBC(ciao, i, mu);
    }
}


void steadyNS::forcesMatrices(label NUmodes, label NPmodes, label NSUPmodes)
{
    tauMatrix.resize(L_U_SUPmodes.size(), 3);
    nMatrix.resize(NPmodes, 3);
    tauMatrix = tauMatrix * 0;
    nMatrix = nMatrix * 0;
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word pName(FORCESdict.lookup("pName"));
    word UName(FORCESdict.lookup("UName"));
    functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);

    for (label i = 0; i < L_U_SUPmodes.size(); i++)
    {
        U = L_U_SUPmodes[i];
        p = Pmodes[0];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            tauMatrix(i, j) = f.forceTau()[j];
        }
    }

    for (label i = 0; i < NPmodes; i++)
    {
        U = L_U_SUPmodes[0];
        p = Pmodes[i];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            nMatrix(i, j) = f.forcePressure()[j];
        }
    }

    if (Pstream::parRun())
    {
        reduce(tauMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::parRun())
    {
        reduce(nMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "eigen", "./ITHACAoutput/Matrices/");
    }
}

void steadyNS::forcesMatrices(label nModes)
{
    tauMatrix.resize(nModes, 3);
    nMatrix.resize(nModes, 3);
    tauMatrix = tauMatrix * 0;
    nMatrix = nMatrix * 0;
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word pName(FORCESdict.lookup("pName"));
    word UName(FORCESdict.lookup("UName"));
    functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);

    for (label i = 0; i < nModes; i++)
    {
        U = Umodes[i];
        p = Pmodes[0];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            tauMatrix(i, j) = f.forceTau()[j];
        }
    }

    for (label i = 0; i < nModes; i++)
    {
        U = Umodes[0];
        p = Pmodes[i];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            nMatrix(i, j) = f.forcePressure()[j];
        }
    }

    if (Pstream::parRun())
    {
        reduce(tauMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::parRun())
    {
        reduce(nMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "eigen", "./ITHACAoutput/Matrices/");
    }
}



void steadyNS::restart()
{
    volScalarField& p = _p();
    volScalarField& p0 = _p0();
    volVectorField& U = _U();
    volVectorField& U0 = _U0();
    surfaceScalarField& phi = _phi();
    surfaceScalarField& phi0 = _phi0();
    p = p0;
    U = U0;
    phi = phi0;
    turbulence.reset(
        (incompressible::turbulenceModel::New(U, phi, _laminarTransport())).ptr()
    );
}

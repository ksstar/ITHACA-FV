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

\*---------------------------------------------------------------------------*/


/// \file
/// Source file of the SteadyBBTurb class.

#include "SteadyBBTurb.H"
#include <cmath>
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
SteadyBBTurb::SteadyBBTurb() {}
SteadyBBTurb::SteadyBBTurb(int argc, char* argv[])
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
#include "createFields.H"
#include "createFvOptions.H"
    //
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "initContinuityErrs.H"
#pragma GCC diagnostic pop
    turbulence->validate();
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
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    viscCoeff = ITHACAdict->lookupOrDefault<word>("viscCoeff", "RBF");
    para = new ITHACAparameters;
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
#include "fvCFD.H"

// Method to performa a truthSolve
void SteadyBB::truthSolve(List<scalar> mu_now, fileName folder)
{ 
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi(); 
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    volScalarField& gh = _gh();
    //volScalarField nut = _nut();
    surfaceScalarField& ghf = _ghf();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
#include "NLsolve.H"
    volScalarField _nut(turbulence->nut());
    //volScalarField nut = turbulence->nut();
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(p_rgh, name(counter), folder);
    ITHACAstream::exportSolution(T, name(counter), folder);
    ITHACAstream::exportSolution(_nut, name(counter), folder);
    Ufield.append(U);
    Pfield.append(p);
    Prghfield.append(p_rgh);
    Tfield.append(T);
    //nutFields.append(_nut);
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
                                   folder);
    }
}

// * * * * * * * * * * * * * * Matrices * * * * * * * * * * * * * * //

List < Eigen::MatrixXd > SteadyBBTurb::turbulenceTerm1(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > ct1Matrix;
    ct1Matrix.setSize(cSize);

    for (label j = 0; j < cSize; j++)
    {
        ct1Matrix[j].resize(nNutModes, cSize);
        ct1Matrix[j] = ct1Matrix[j] * 0;
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct1Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                         nutModes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(ct1Matrix, "ct1Matrix", "eigen",
                               "./ITHACAoutput/Matrices/ct1");
    return ct1Matrix;
}

Eigen::Tensor<double, 3> SteadyBBTurb::turbulenceTensor1(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor;
    ct1Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                         nutModes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct1Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
                                  "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct1Tensor;
}

List < Eigen::MatrixXd > SteadyBBTurb::turbulenceTerm2(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > ct2Matrix;
    ct2Matrix.setSize(cSize);

    for (label j = 0; j < cSize; j++)
    {
        ct2Matrix[j].resize(nNutModes, cSize);
        ct2Matrix[j] = ct2Matrix[j] * 0;
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct2Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(ct2Matrix, "ct2Matrix", "eigen",
                               "./ITHACAoutput/Matrices/ct2");
    return ct2Matrix;
}

Eigen::Tensor<double, 3> SteadyBBTurb::turbulenceTensor2(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor;
    ct2Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                  "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

Eigen::MatrixXd SteadyBBTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    label btSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
    btMatrix = btMatrix * 0;

    // Project everything
    for (label i = 0; i < btSize; i++)
    {
        for (label j = 0; j < btSize; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(dev((T(
                    fvc::grad(
                        L_U_SUPmodes[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
                                  "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return btMatrix;
}



// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void SteadyBBTurb::projectSUP(fileName folder, label NU, label NP, label NT,
                            label NSUP, label Nnut)
{
    NUmodes = NU;
    NPrghmodes = NP;
    NTmodes = NT;
    NSUPmodes = NSUP;
    nNutModes = Nnut;

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        
        word bStr = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bStr))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", bStr);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
        {
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        }
        else
        {
            btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        word cStr = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + cStr))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", cStr);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPrghmodes, NSUPmodes);
        }

	word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        }
        else
        {
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        }
        else
        {
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        word kStr = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + kStr))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", kStr);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word hStr = "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) +  "_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + hStr))
        {
            ITHACAstream::ReadDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/", hStr);
        }
        else
        {
            H_matrix = buoyant_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word wStr = "W_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + wStr))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", wStr);
        }
        else
        {
            W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

	word qStr = "Q_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + qStr))
        {
            ITHACAstream::ReadDenseMatrix(Q_matrix, "./ITHACAoutput/Matrices/", qStr);
        }
        else
        {
	    Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);            
        }

        word yStr = "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + yStr))
        {
            ITHACAstream::ReadDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/", yStr);
        }
        else
        {
            Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

        word pStr = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + pStr))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", pStr);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
        }


	/*if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }*/
    }
    else
    {
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

        L_T_modes.resize(0);

        if (liftfieldT.size() != 0)
        {
            for (label k = 0; k < liftfieldT.size(); k++)
            {
                L_T_modes.append(liftfieldT[k]);
            }
        }

        if (NTmodes != 0)
        {
            for (label k = 0; k < NTmodes; k++)
            {
                L_T_modes.append(Tmodes[k]);
            }
        }

        M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        H_matrix = buoyant_term(NUmodes, NTmodes, NSUPmodes);
        W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);

	/*if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }*/
    }


    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields,
              nutModes, nNutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                               "./ITHACAoutput/Matrices/");
    // Export the matrix
    ITHACAstream::SaveDenseMatrix(coeffL2, "./ITHACAoutput/Matrices/",
                                  "coeffL2_nut_" + name(nNutModes));
    samples.resize(nNutModes);
    rbfSplines.resize(nNutModes);
    Eigen::MatrixXd weights;

    for (int i = 0; i < nNutModes; i++)
    {
        word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                          + name(NUmodes) + "_" + name(NSUPmodes) ;

        if (ITHACAutilities::check_file("./ITHACAoutput/weightsSUP/" + weightName))
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/",
                                          weightName);
            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
        else
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN);
            ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                          "./ITHACAoutput/weightsSUP/", weightName);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
    }
}



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
void SteadyBBTurb::truthSolve(List<scalar> mu_now, fileName folder)
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
    surfaceScalarField& ghf = _ghf();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    volScalarField& nut = _nut();
    volScalarField& k = _k();
#include "NLsolve.H"
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(p_rgh, name(counter), folder);
    ITHACAstream::exportSolution(T, name(counter), folder);
    ITHACAstream::exportSolution(nut, name(counter), folder);
    ITHACAstream::exportSolution(k, name(counter), folder);
    Ufield.append(U);
    Pfield.append(p);
    Prghfield.append(p_rgh);
    Tfield.append(T);
    nutFields.append(nut);
    kFields.append(k);
    counter++;
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
                    fvc::grad(L_U_SUPmodes[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
                                  "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));

    ITHACAstream::exportMatrix(btMatrix, "btMatrix", "eigen", "./ITHACAoutput/Matrices/");

    return btMatrix;
}

List <Eigen::MatrixXd> SteadyBBTurb::temperature_turbulence_term(
    label NTmodes, label nNutModes)
{
    label Stsize = NTmodes + liftfieldT.size();
    List <Eigen::MatrixXd> S_matrix;
    S_matrix.setSize(Stsize);

    for (label j = 0; j < Stsize; j++)
    {
        S_matrix[j].resize(nNutModes, Stsize);
    }

    for (label i = 0; i < Stsize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix S_matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < Stsize; k++)
            {
                S_matrix[i](j, k) = fvc::domainIntegrate(L_T_modes[i] * fvc::laplacian(
                                        nutModes[j], L_T_modes[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "eigen",
                               "./ITHACAoutput/Matrices/S");
    return S_matrix;
}

/*
Eigen::MatrixXd SteadyBBTurb::pressure_turbulence_term(label NUmodes, label nTkeModes)
{
    label pt1size = NUmodes + NSUPmodes + liftfield.size();
    label pt2size = nTkeModes;
    Eigen::MatrixXd ptMatrix(pt1size, pt2size);

    // Project everything
    for (label i = 0; i < pt1size; i++)
    {
        for (label j = 0; j < pt2size; j++)
        {
            ptMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::grad(
                    kModes[j])).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(ptMatrix, "./ITHACAoutput/Matrices/",
                                  "pt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nTkeModes));
    return ptMatrix;
}


Eigen::MatrixXd SteadyBBTurb::pressure_term(label NPmodes, label NPrghmodes)
{
    label pt1size = NPrghmodes ;
    label pt2size = NPmodes;
    Eigen::MatrixXd rMatrix(pt1size, pt2size);

    // Project everything
    for (label i = 0; i < pt1size; i++)
    {
        for (label j = 0; j < pt2size; j++)
        {
            rMatrix(i, j) = fvc::domainIntegrate(Prghmodes[i] *Pmodes[j]).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(rMatrix, "./ITHACAoutput/Matrices/",
                                  "R_" + name(NPmodes) + "_" + name(NPrghmodes));
    return rMatrix;
}

Eigen::MatrixXd SteadyBBTurb::pressure_term2(label NPmodes, label NPrghmodes)
{
    label pt1size = NPrghmodes ;
    label pt2size = NPrghmodes;
    Eigen::MatrixXd r2Matrix(pt1size, pt2size);

    // Project everything
    for (label i = 0; i < pt1size; i++)
    {
        for (label j = 0; j < pt2size; j++)
        {
            r2Matrix(i, j) = fvc::domainIntegrate(Prghmodes[i] *Prghmodes[j]).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(r2Matrix, "./ITHACAoutput/Matrices/",
                                  "R2_" + name(NPmodes) + "_" + name(NPrghmodes));
    return r2Matrix;
}

Eigen::MatrixXd SteadyBBTurb::pressure_term3(label NPmodes, label NPrghmodes)
{
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedVector& g = _g();
    surfaceScalarField& ghf = _ghf();
    volScalarField& gh = _gh();

    label pt1size = NPrghmodes ;
    label pt2size = NTmodes;
    Eigen::MatrixXd r3Matrix(pt1size, pt2size);

    // Project everything
    for (label i = 0; i < pt1size; i++)
    {
        for (label j = 0; j < pt2size; j++)
        {
            r3Matrix(i, j) = fvc::domainIntegrate(Prghmodes[i] *(gh*(1.0 - (beta * (L_T_modes[j] - TRef))))).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(r3Matrix, "./ITHACAoutput/Matrices/",
                                  "R3_" + name(NPmodes) + "_" + name(NPrghmodes));
    return r3Matrix;
}*/

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //
Eigen::MatrixXd SteadyBBTurb::pressure_gradient_term2(label NUmodes,
        label NPmodes,
        label NSUPmodes)
{
    label K1size = NUmodes + NSUPmodes + liftfield.size();
    label K2size = NPmodes;
    Eigen::MatrixXd K2_matrix(K1size, K2size);

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K2_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                                                  fvc::reconstruct(fvc::snGrad(Pmodes[j]) *
                                                          Pmodes[j].mesh().magSf())).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(K2_matrix, "./ITHACAoutput/Matrices/",
                                  "K2_" + name(liftfield.size()) + "_" + name(NUmodes)
                                  + "_" + name(NSUPmodes) + "_" + name(NPmodes));
    return K2_matrix;
}

/*
Eigen::MatrixXd SteadyBBTurb::pressure_gradient_term2(label NPmodes,
        label NPrghmodes)
{


    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix2(K2size, K2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < K2size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Prghmodes[i]);
	    surfaceScalarField BC4 = fvc::snGrad(Pmodes[i]);
	    surfaceScalarField BC5 = (BC3 * BC4) * mesh.magSf();

            double s = 0;

	    for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                    s += gSum(BC5.boundaryField()[k]);
            }

            K_matrix2(i, j) = s;
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(K_matrix2, "./ITHACAoutput/Matrices/",
                                  "K2_" + name(NPmodes)
                                  + "_" + name(NPmodes));
    return K_matrix2;
}

Eigen::MatrixXd SteadyBBTurb::pressure_gradient_term3(label NPmodes,
        label NPrghmodes)
{
    label K1size = NPmodes;
    label K2size = NPrghmodes;
    Eigen::MatrixXd K_matrix3(K1size, K2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Prghmodes[i]);
	    surfaceScalarField BC4 = fvc::snGrad(Prghmodes[j]);
	    surfaceScalarField BC5 = (BC3 * BC4) * mesh.magSf();

            double s = 0;

	    for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                    s += gSum(BC5.boundaryField()[k]);
            }

            K_matrix3(i, j) = s;
        }
    }
 
    return K_matrix3;
}

Eigen::MatrixXd SteadyBBTurb::pressure_gradient_term4(label NPmodes,
        label NPrghmodes)
{


    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix4(K2size, K2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < K2size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Pmodes[i]);
	    surfaceScalarField BC4 = fvc::snGrad(Pmodes[i]);
	    surfaceScalarField BC5 = (BC3 * BC4) * mesh.magSf();

            double s = 0;

	    for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                    s += gSum(BC5.boundaryField()[k]);
            }

            K_matrix4(i, j) = s;
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(K_matrix4, "./ITHACAoutput/Matrices/",
                                  "K4_" + name(NPmodes)
                                  + "_" + name(NPmodes));
    return K_matrix4;
}


Eigen::MatrixXd SteadyBBTurb::buoyant_term2(label NPmodes, label NPrghmodes, label NTmodes)
{
    label H1size = NPrghmodes;
    label H2size = NTmodes + liftfieldT.size() ;
    Eigen::MatrixXd H_matrix2(H1size, H2size);
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedVector& g = _g();
    surfaceScalarField& ghf = _ghf();
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Prghmodes[i]);
	    surfaceScalarField BC4 = fvc::snGrad(1.0 - (beta * (L_T_modes[j] - TRef)));
	    surfaceScalarField BC5 = (BC3 * (ghf * BC4)) * mesh.magSf();

            double s = 0;

	    for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                    s += gSum(BC5.boundaryField()[k]);
            }

            H_matrix2(i, j) = s;
        }
    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(H_matrix2, "./ITHACAoutput/Matrices/",
                                  "H2_" + name(NPmodes) +  "_" + name(NTmodes));
    return H_matrix2;
}



Eigen::MatrixXd SteadyBBTurb::buoyant_term3(label NPmodes, label NPrghmodes, label NTmodes)
{
    label H1size = NPrghmodes;
    label H2size = NTmodes + liftfieldT.size() ;
    Eigen::MatrixXd H_matrix3(H1size, H2size);
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedVector& g = _g();
    surfaceScalarField& ghf = _ghf();
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
	    H_matrix3(i, j) = fvc::domainIntegrate(Pmodes[i] *(9.81*(- (beta * (L_T_modes[j] - TRef))))).value();
        }
    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(H_matrix3, "./ITHACAoutput/Matrices/",
                                  "H3_" + name(NPmodes) +  "_" + name(NTmodes));
    return H_matrix3;
}*/


Eigen::MatrixXd SteadyBBTurb::buoyant_term4(label NUmodes, label NTmodes, label NSUPmodes)
{
    label H1size = NUmodes + NSUPmodes + liftfield.size();;
    label H2size = NTmodes + liftfieldT.size() ;
    Eigen::MatrixXd H_matrix4(H1size, H2size);
  /*  dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedVector& g = _g();
    surfaceScalarField& ghf = _ghf();
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
	  //  H_matrix4(i, j) = fvc::domainIntegrate(fvc::reconstruct(
		//			((n & fvc::interpolate(L_U_SUPmodes[i]))
		//		))).value();


            //H_matrix4(i, j) = fvc::domainIntegrate(fvc::reconstruct( fvc::interpolate(Tmodes[i]) * Tmodes[i].mesh().magSf())).value()
//* 9.81*beta *fvc::interpolate( L_T_modes[j])

             H_matrix4(i, j)  = fvc::domainIntegrate(L_U_SUPmodes[i] & (g*fvc::reconstruct(fvc::interpolate(Tmodes[j]) *
                                                          Tmodes[j].mesh().magSf()))).value();

        }
    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(H_matrix4, "./ITHACAoutput/Matrices/",
                                  "H4_" + name(NUmodes) +  "_" + name(NTmodes));*/
    return H_matrix4;
}




Eigen::MatrixXd SteadyBBTurb::buoyant_term5(label NUmodes, label NTmodes, label NSUPmodes)
{
    label H1size = NUmodes + NSUPmodes + liftfield.size();;
    Eigen::MatrixXd H_matrix5(H1size, 1);
    /*dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedVector& g = _g();
    surfaceScalarField& ghf = _ghf();
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < H1size; i++)
    {

             H_matrix5(i,0)  =  fvc::domainIntegrate(L_U_SUPmodes[i] & g).value();

    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(H_matrix5, "./ITHACAoutput/Matrices/",
                                  "H5_" + name(NUmodes) +  "_" + name(NTmodes));*/
    return H_matrix5;
}


// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void SteadyBBTurb::projectSUP(fileName folder, label NU, label NP, label NPrgh, label NT,
                            label NSUP, label Nnut, label NTKE, Eigen::MatrixXd para2)
{
    NUmodes = NU;
    NPrghmodes = NPrgh;
    NPmodes = NP;
    NTmodes = NT;
    NSUPmodes = NSUP;
    nNutModes = Nnut;
    nTkeModes = NTKE;

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        
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

        B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
	//K2_matrix = pressure_gradient_term2(NUmodes, NPmodes, NSUPmodes);
        H_matrix = buoyant_term(NUmodes, NTmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
 	btMatrix = btTurbulence(NUmodes, NSUPmodes);
	ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
	ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
	S_matrix = temperature_turbulence_term(NTmodes, nNutModes);
	//H_matrix4 = buoyant_term4(NUmodes, NTmodes, NSUPmodes);
	//H_matrix5 = buoyant_term5(NUmodes, NTmodes, NSUPmodes);

	if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
	    bcTempVec = bcTemperatureVec(NTmodes);
            bcTempMat = bcTemperatureMat(NTmodes);
	    bcGradTVec = bcGradientTVec(NTmodes);
            bcGradTMat = bcGradientTMat(NTmodes);
        }
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes;
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
            	samples[i]->addSample(para2.row(j), coeffL2(i, j)); // samples[i]->addSample(mu.row(j), coeffL2(i, j));
           }

	ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/",
                                          weightName);

        rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, e);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    	}
    	else
    	{
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(para2.row(j), coeffL2(i, j));
            }

            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, 0, e);
            ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                          "./ITHACAoutput/weightsSUP/", weightName);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
    }



    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
  /*  tkeCoeffL2 = ITHACAutilities::get_coeffs_ortho(kFields,
              kModes, nTkeModes);

    tkeSamples.resize(nTkeModes);
    tkeRbfSplines.resize(nTkeModes);
    
     for (int i = 0; i < nTkeModes; i++)
    {
        tkeSamples[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < tkeCoeffL2.cols(); j++)
        {
	     Info << "tkeSamples"  << endl;
            tkeSamples[i]->addSample(para.row(j), tkeCoeffL2(i, j)); // samples[i]->addSample(mu.row(j), coeffL2(i, j));
        }

        tkeRbfSplines[i] = new SPLINTER::RBFSpline(*tkeSamples[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction TKE for mode " << i + 1 << std::endl;
    }*/



}



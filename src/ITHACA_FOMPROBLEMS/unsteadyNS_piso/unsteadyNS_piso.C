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
/// Source file of the unsteadyNS class.

#include "unsteadyNS_piso.H"
#include "viscosityModel.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
unsteadyNS_piso::unsteadyNS_piso() {}

unsteadyNS_piso::unsteadyNS_piso(int argc, char* argv[])
:
unsteadyNS_PISO(argc, argv)
{
  Info << offline << endl; 
}

fvVectorMatrix unsteadyNS_piso::get_Umatrix(volVectorField& U,
  volScalarField& p)
{
#include "initContinuityErrs.H"
  Time& runTime = _runTime();
  fvMesh& mesh = _mesh();
  IOMRFZoneList& MRF = _MRF();
  surfaceScalarField& phi = _phi();
  fv::options& fvOptions = _fvOptions();
  MRF.correctBoundaryVelocity(U);

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

  dimensionedScalar nu
  (
   transportProperties.lookup("nu")
   );

  fvVectorMatrix Ueqn
  (
    fvm::ddt(U) 
    + fvm::div(phi, U)
    -fvm::laplacian(nu,U)
    );
  Ueqn_global = &Ueqn;
  return Ueqn;
}

fvScalarMatrix unsteadyNS_piso::get_Pmatrix(volVectorField& U,
  volScalarField& p)
{
  Time& runTime = _runTime();
  IOMRFZoneList& MRF = _MRF();
  surfaceScalarField& phi = _phi();
  MRF.correctBoundaryVelocity(U);
  fvMesh& mesh = _mesh();
  pisoControl piso(mesh);

  volScalarField rAU(1.0 / Ueqn_global->A());
  volVectorField HbyA(constrainHbyA(rAU * Ueqn_global->H(), U, p));
  surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA)+  		  fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));
  adjustPhi(phiHbyA, U, p);

    // Update the pressure BCs to ensure flux consistency    
  constrainPressure(p, U, phiHbyA, rAU); 

  while (piso.correctNonOrthogonal()) 
  { 
    fvScalarMatrix pEqn
    (
      fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
      );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter()))); 
    
    if (piso.finalNonOrthogonalIter())
    {
      phi = phiHbyA - pEqn.flux();
    }
  }
    #include "continuityErrs.H"
  U = HbyA - rAU * fvc::grad(p);
  U.correctBoundaryConditions();

  fvScalarMatrix pEqn
  (
    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );
  return pEqn;
}

void unsteadyNS_piso::truthSolve2(List<scalar> mu_now)
{
  Time& runTime = _runTime();
  IOMRFZoneList& MRF = _MRF();
  volScalarField& p = _p();
  volVectorField& U = _U();
  surfaceScalarField& phi = _phi();
  fvMesh& mesh = _mesh();
  fv::options& fvOptions = _fvOptions();
  pisoControl piso(mesh);

  singlePhaseTransportModel& laminarTransport = _laminarTransport();
  instantList Times = runTime.times();
  runTime.setEndTime(finalTime);
    // Perform a TruthSolve
  runTime.setTime(Times[1], 1);
  runTime.setDeltaT(timeStep);
  nextWrite = startTime;
    // save initial condition in folder 0
  ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
  ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
  std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
   runTime.timeName());
  Ufield.append(U);
  Pfield.append(p);
  counter++;
    //nextWrite += writeEvery;


   // Start the time loop
  while (runTime.loop())
  {

   Info << "Time = " << runTime.timeName() << nl << endl;
	#include "CourantNo.H"

   fvVectorMatrix UEqn(get_Umatrix(U, p));

   if (piso.momentumPredictor())
   {
     solve(UEqn == -fvc::grad(p));
   }

            // --- Pressure corrector loop
   while (piso.correct())
   {
    get_Pmatrix(U, p);
  }

  Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
  << nl << endl;

       // if (checkWrite(runTime))
       // {
  ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
  ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    //        std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +  runTime.timeName());
  Ufield.append(U);
  Pfield.append(p);
  counter++;

}

}

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

#include "ReducedSteadyBBTurb.H"

#include "fvc.H"
//#include "fvm.H"
//#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedSteadyBBTurb::ReducedSteadyBBTurb()
{
}

ReducedSteadyBBTurb::ReducedSteadyBBTurb(SteadyBBTurb& fomProblem)
    :
    problem(&fomProblem)
{
    N_BC_t    = problem->inletIndexT.rows();
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = 0;
    Nphi_prgh = 0;
    //Nphi_prgh = problem->K_matrix.cols();
    Nphi_t    = problem->Y_matrix.rows();	
    nphiNut = problem->cTotalTensor.dimension(1);

    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        LUmodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        LUmodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        LUmodes.append(problem->supmodes[k]);
    }

    for (label k = 0; k < problem->NTmodes; k++)
    {
        LTmodes.append(problem->Tmodes[k]);
    }

  /*  for (label k = 0; k < problem->NPmodes; k++)
    {
        LPmodes.append(problem->Pmodes[k]);
    }

   for (label k = 0; k < problem->NPrghmodes; k++)
    {
        LPrghmodes.append(problem->Prghmodes[k]);
    }*/

    newtonObject = newtonSteadyBBTurb(Nphi_u + Nphi_prgh + Nphi_t, Nphi_u + Nphi_prgh + Nphi_t, fomProblem);
}

int newtonSteadyBBTurb::operator()(const Eigen::VectorXd& x,
                                   Eigen::VectorXd& fvec) const
{
//std::cout << "HERE 1 "  << std::endl;
    Eigen::VectorXd aTmp(Nphi_u);
    //Eigen::VectorXd bTmp(Nphi_prgh);
    Eigen::VectorXd cTmp(Nphi_t);
    aTmp = x.head(Nphi_u);
   // bTmp = x.segment(Nphi_u, Nphi_prgh);
    cTmp = x.tail(Nphi_t);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd dd(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu ;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * aTmp;
    //Eigen::VectorXd m2 = problem->K_matrix * pCnst;
    // Pressure Term
    //Eigen::VectorXd m3 = problem->P_matrix * aTmp;
    // Buoyancy Term
    Eigen::VectorXd m10 = problem->H_matrix * cTmp;
    //Eigen::VectorXd m11 = problem->H_matrix4 * cTmp*2.5644e-04;
    //Eigen::VectorXd m12 = problem->H_matrix5.col(0)*(1-2.5644e-04*423.15);
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    Eigen::MatrixXd st(1, 1);
    Eigen::MatrixXd st2(1, 1);
    // diffusive term temperature
    Eigen::VectorXd m6 = problem->Y_matrix * cTmp * (nu / Pr) ;
    // Penalty term - Temperature
    Eigen::MatrixXd penaltyT = Eigen::MatrixXd::Zero(Nphi_t, N_BC_t+1);
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    //Vector for penalty method to enforce spatial dependent velocity boundary condition
    Eigen::MatrixXd bcVelVec2 = Eigen::MatrixXd::Zero(Nphi_u, 1);
    Eigen::VectorXd Out;
    for (label l = 0; l < problem->NUmodes; l++)
    {
    	unsigned int size = problem->Umodes[l].boundaryField()[0].size();
        Out.resize(size);
	for (unsigned int k = 0; k < size ; k++)
        {
             Out(k) =  problem->Umodes[l].boundaryField()[0][k][0];
        }
	 bcVelVec2(l, 0) = bc.dot(Out);
    }

   /* for (label l = 0; l < problem->NSUPmodes; l++)
    {
	label j = l + problem->NUmodes;
    	unsigned int size = problem->supmodes[l].boundaryField()[0].size();
        Out.resize(size);
	for (unsigned int k = 0; k < size ; k++)
        {
             Out(k) =  problem->supmodes[l].boundaryField()[0][k][0];
        }
	bcVelVec2(j, 0) = bc.dot(Out);
    }*/

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (label l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l) * (bcVelVec2 - problem->bcVelMat[l] *
                              aTmp);
	}

	for (label l = 0; l < N_BC_t; l++)
        {
 	
	    penaltyT.col(l) = tauT(l) * (bc_t(l) *  problem->bcTempVec[l] - problem->bcTempMat[l] * cTmp);
	   
        }

	penaltyT.col(N_BC_t) =  tauHF(0) *( gradient_t(0,0) * problem->bcGradTVec[0] - problem->bcGradTMat[0] * cTmp); 
    }

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp -gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp ;
        fvec(i) = m1(i) - cc(0, 0)  - m2(i)-m10(i);
	if (problem->bcMethod == "penalty")
        {
            for (label l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }

    }

   // for (label j = 0; j < Nphi_prgh; j++)
   // {
        //label k = j + Nphi_u;
	//fvec(k) = m3(j);
	//fvec(k) = pCnst(j);
   // }
//std::cout << "HERE 4 "  << std::endl;
    for (label i = 0; i < Nphi_t; i++)
    {
	//label k = i + Nphi_u + Nphi_prgh;
	label k = i + Nphi_u ;
        qq = aTmp.transpose() * problem->Q_matrix[i] * cTmp;
        //st = gNut.transpose() * problem->S_matrix[i] * cTmp;
	st2 = gAlphat.transpose() * problem->S2_matrix[i] * cTmp;

        //fvec(k) = m6(i) - qq(0, 0) + st(0, 0) / Prt;
	fvec(k) = m6(i) - qq(0, 0) + st2(0, 0) ;

 	if (problem->bcMethod == "penalty")
        {
            for (label l = 0; l < N_BC_t+1; l++)
            {
                fvec(k) += penaltyT(i, l);
            }
        }

    }	
//td::cout << "HERE 5 "  << std::endl;
    if (problem->bcMethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - bc(j);
        }

	for (label j = 0; j < N_BC_t; j++)
    	{
	    label k = j + Nphi_u + Nphi_prgh;
            fvec(j) = x(k) - bc_t(j);
    	}
    }

//std::cout << "fvec = " << fvec << std::endl;

    return 0;
}

int newtonSteadyBBTurb::df(const Eigen::VectorXd& x,
                           Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonSteadyBBTurb> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //


void ReducedSteadyBBTurb::solveOnlineSUP(Eigen::MatrixXd vel, Eigen::MatrixXd temp, int param )
{

    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();

    //scalar param2 = param;
    //scalar index;
    Eigen::MatrixXf::Index index;
    // find nearest neighbour
    (parOn.array()-gradient_t(0,0)).cwiseAbs().minCoeff(&index);

    std::cout << "index = " << index << std::endl;

     y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[index],
                     LUmodes);

    y.tail(Nphi_t) = ITHACAutilities::get_coeffs(problem->Tfield[index],
                     LTmodes);

    if (problem->bcMethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel(j, 0);
        }

        for (label j = 0; j < N_BC_t; j++)
        {
            label k = j + Nphi_prgh + Nphi_u;
            y(k) = temp(j, 0);
        } 
    }



    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newtonSteadyBBTurb> hnls(newtonObject);
    newtonObject.bc = vel.col(0);
    newtonObject.bc_t.resize(N_BC_t);
    newtonObject.tauU = tauU;
    newtonObject.tauT = tauT;
    newtonObject.tauHF = tauHF; 
    newtonObject.gradient_t = gradient_t;
    newtonObject.nu = nu;
    newtonObject.Pr = Pr;
    newtonObject.Prt = Prt;
    //newtonObject.pCnst = ITHACAutilities::get_coeffs(problem->Prghfield[param],
      //               LPrghmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC_t; j++)
    {
        newtonObject.bc_t(j) = temp(j, 0);
    }

rbfCoeff.resize(nphiNut);

    // Turbulence 
//    if (problem->viscCoeff == "L2")
 //   {
 //       for (label i = 0; i < nphiNut; i++)
 //       {
 //           newtonObject.gNut = problem->nutCoeff;
 //       }
 //   }
 //   else if (problem->viscCoeff == "RBF")
 //   {

std::cout << "gradient_t = " << gradient_t << std::endl;

        for (label i = 0; i < nphiNut; i++)
        {
            newtonObject.gNut(i) = problem->rbfSplines[i]->eval(gradient_t);
	    rbfCoeff(i) = newtonObject.gNut(i);
        }

std::cout << "newtonObject.gNut(0) = " << newtonObject.gNut(0) << std::endl;

       // rbfCoeffMat(0, param) = param;
       // rbfCoeffMat.block(1, param, nphiNut, 1) = newtonObject.gNut;

 //   }
//    else
//    {
//        Info << "The way to compute the eddy viscosity coefficients has to be either L2 or RBF"
//             << endl;
//        exit(0);
//    }

    volScalarField nutTmp("nutRec", problem->nutModes[0] * 0);

    for (label j = 0; j < nphiNut; j++)
    {
        nutTmp += problem->nutModes[j] * newtonObject.gNut(j);
    }
    ITHACAstream::exportSolution(nutTmp, name(count_online_solve), "./ITHACAoutput/ReconstructionSUP");
    nutRec.append(nutTmp);

   for (label i = 0; i < nphiNut; i++)
        {
            newtonObject.gAlphat(i) = problem->alphatRbfSplines[i]->eval(gradient_t);
	    rbfCoeff(i) = newtonObject.gAlphat(i);
        }

    volScalarField alphatTmp("alphatRec", problem->alphatModes[0] * 0);

    for (label j = 0; j < nphiNut; j++)
    {
        alphatTmp += problem->alphatModes[j] * newtonObject.gAlphat(j);
    }
    ITHACAstream::exportSolution(alphatTmp, name(count_online_solve), "./ITHACAoutput/ReconstructionSUP");
    alphatRec.append(alphatTmp);

    hnls.solve(y);
//std::cout << "HERE C "  << std::endl;
    Eigen::VectorXd res(y);
    newtonObject.operator()(y, res);
    std::cout << "################## Online solve N° " << count_online_solve <<
              " ##################" << std::endl;

    if (res.norm() < 1e-5)
    {
        std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    else
    {
        std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }

    count_online_solve += 1;

//exit(0);
}


void ReducedSteadyBBTurb::reconstructSUP(fileName folder, int printEvery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            volVectorField uRec("uRec", LUmodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                uRec += LUmodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(uRec, name(online_solution[i](0, 0)), folder);

            volScalarField prghRec("prghRec", problem->Prghmodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                prghRec += problem->Prghmodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(prghRec, name(online_solution[i](0, 0)), folder);


	    volScalarField TRec("TRec", problem->Tmodes[0] * 0);

            for (label j = 0; j < Nphi_t; j++)
            {
                TRec += problem->Tmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

	    ITHACAstream::exportSolution(TRec, name(online_solution[i](0, 0)), folder);

	    scalar patchHeatFlux2 = gSum((TRec.boundaryFieldRef()[4] - TRec.boundaryFieldRef()[4].patchInternalField())*
			(TRec.mesh().deltaCoeffs().boundaryField()[4])) ;

	    std::cout << "patchHeatFlux2: " << patchHeatFlux2 << std::endl;

	   // scalar patchHeatFlux2F = gSum((problem->Tfield[i].boundaryFieldRef()[4] - problem->Tfield[i].boundaryFieldRef()[4].patchInternalField())*
			//(problem->Tfield[i].mesh().deltaCoeffs().boundaryField()[4])) ;

	   //std::cout << "patchHeatFluxFOM: " << patchHeatFlux2F << std::endl;

            nextWrite += printEvery;
            UREC.append(uRec);
	    TREC.append(TRec);
            PREC.append(prghRec);

	// Calculation local Stanton number
    	   label stantonSize = TRec.boundaryFieldRef()[4].size();

	   //Eigen::MatrixXd stantonOLD(stantonSize, 1);
    	   Eigen::MatrixXd stanton(stantonSize, 1);
	   Eigen::MatrixXd CoordHeater(stantonSize, 1);
	   Eigen::MatrixXd Cf(stantonSize, 1);

	   vectorField BC5 = uRec.boundaryField()[4].snGrad();
	   for (label j = 0; j < stantonSize; j++)
	   {
	   	//stantonOLD(j,0) = (gradient_t(0,0) * h) / (lambda * (TRec.boundaryField()[4][j] - Tref));
		//stantonOLD(j,0) = stantonOLD(j,0)/ (Re * Pr);
		//stanton(j,0) = (gradient_t(0,0)*lambda*Pr) / (Ub *nu * (TRec.boundaryField()[4][j] - Tref));
		stanton(j,0) = (gradient_t(i,0)*nu) / (Ub *lambda*Pr* (TRec.boundaryField()[4][j] - Tref));
		CoordHeater(j,0) = TRec.mesh().boundaryMesh()[4].faceCentres()[j].x();
		Cf(j,0) = (2*nu*-BC5[j].component(0))/(Ub*Ub)*1000;

	   }

	   ITHACAstream::exportMatrix(stanton, "stantonROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   //ITHACAstream::exportMatrix(stantonOLD, "stantonOLDROM", "eigen",
                         //      "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(CoordHeater, "CoordHeaterROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(Cf, "CfROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));


	   label stantonLWSize = TRec.boundaryFieldRef()[3].size();
	   Eigen::MatrixXd stantonLW(stantonLWSize, 1);
 	  // Eigen::MatrixXd stantonLWOLD(stantonLWSize, 1);
	   Eigen::MatrixXd CoordLW(stantonSize, 1);
	   Eigen::MatrixXd CfLW(stantonSize, 1);
	   vectorField BC5LW = uRec.boundaryField()[3].snGrad();

	   for (label j = 0; j < stantonLWSize; j++)
	   {
	   //	stantonLWOLD(j,0) = (gradient_t(0,0) * h) / (lambda * (TRec.boundaryField()[3][j] - Tref));
		//stantonLWOLD(j,0) = stantonLWOLD(j,0)/ (Re * Pr);
		//stantonLW(j,0) = (gradient_t(0,0)*lambda*Pr) / (Ub *nu * (TRec.boundaryField()[3][j] - Tref));
		stantonLW(j,0) = (gradient_t(i,0)*nu) / (Ub *lambda*Pr* (TRec.boundaryField()[3][j] - Tref));
		CoordLW(j,0) = TRec.mesh().boundaryMesh()[3].faceCentres()[j].x();
		CfLW(j,0) = (2*nu*-BC5LW[j].component(0))/(Ub*Ub)*1000;
	   }

	  ITHACAstream::exportMatrix(stantonLW, "stantonLWROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	//ITHACAstream::exportMatrix(stantonLWOLD, "stantonLWROMOLD", "eigen",
                  //             "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(CoordLW, "CoordLWROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	ITHACAstream::exportMatrix(CfLW, "CfLWROM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));
    
        }

        counter++;
    }
}


void ReducedSteadyBBTurb::nonDimNum()
{
    for (label i = 0; i < problem->Ufield.size(); i++)
    {
    	   label stantonSize = problem->Tfield[i].boundaryFieldRef()[4].size();

    	   Eigen::MatrixXd stanton(stantonSize, 1);
	   Eigen::MatrixXd stantonOLD(stantonSize, 1);
	   Eigen::MatrixXd CoordHeater(stantonSize, 1);
	   Eigen::MatrixXd Cf(stantonSize, 1);
	   vectorField BC5 = problem->Ufield[i].boundaryField()[4].snGrad();

	   for (label j = 0; j < stantonSize; j++)
	   {
	   	stantonOLD(j,0) = (gradient_t(i,0) * h) / (lambda * (problem->Tfield[i].boundaryField()[4][j] - Tref));
		stantonOLD(j,0) = stantonOLD(j,0)/ (Re * Pr);
	        stanton(j,0) = (gradient_t(i,0)*lambda*Pr) / (Ub *nu * (problem->Tfield[i].boundaryField()[4][j] - Tref));
		CoordHeater(j,0) = problem->Tfield[i].mesh().boundaryMesh()[4].faceCentres()[j].x();
		Cf(j,0) = (2*nu*-BC5[j].component(0))/(Ub*Ub)*1000;

	   }

	   ITHACAstream::exportMatrix(stanton, "stantonFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(stantonOLD, "stantonFOMOLD", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(CoordHeater, "CoordHeaterFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(Cf, "CfFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));


	   label stantonLWSize = problem->Tfield[i].boundaryFieldRef()[3].size();
	   Eigen::MatrixXd stantonLW(stantonLWSize, 1);
	   Eigen::MatrixXd stantonLWOLD(stantonLWSize, 1);
	   Eigen::MatrixXd CoordLW(stantonSize, 1);
	   Eigen::MatrixXd CfLW(stantonSize, 1);
	   vectorField BC5LW = problem->Ufield[i].boundaryField()[3].snGrad();

	   for (label j = 0; j < stantonLWSize; j++)
	   {
	   	stantonLWOLD(j,0) = (gradient_t(i,0) * h) / (lambda * (problem->Tfield[i].boundaryField()[3][j] - Tref));
		stantonLWOLD(j,0) = stantonLWOLD(j,0)/ (Re * Pr);
	        stantonLW(j,0) = (gradient_t(i,0)*lambda*Pr) / (Ub *nu * (problem->Tfield[i].boundaryField()[3][j] - Tref));
		CoordLW(j,0) = problem->Tfield[i].mesh().boundaryMesh()[3].faceCentres()[j].x();
		CfLW(j,0) = (2*nu*-BC5LW[j].component(0))/(Ub*Ub)*1000;
	   }


	  ITHACAstream::exportMatrix(stantonLW, "stantonLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(stantonLWOLD, "stantonLWFOMOLD", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(CoordLW, "CoordLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(CfLW, "CfLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

    }

}


void ReducedSteadyBBTurb::nonDimNumOn()
{


    for (label i = 0; i < problem->Ufield_on.size(); i++)
    {
    	   label stantonSize = problem->Tfield_on[i].boundaryFieldRef()[4].size();

    	   Eigen::MatrixXd stanton(stantonSize, 1);
	   //Eigen::MatrixXd stantonOLD(stantonSize, 1);
	   Eigen::MatrixXd CoordHeater(stantonSize, 1);
	   Eigen::MatrixXd Cf(stantonSize, 1);
	   vectorField BC5 = problem->Ufield_on[i].boundaryField()[4].snGrad();

std::cout << "### gradient(i,0) " << gradient_t(i,0)<< " ###" << std::endl; 

	   for (label j = 0; j < stantonSize; j++)
	   {
		
	   	//stantonOLD(j,0) = (gradient_t(i,0) * h) / (lambda * (problem->Tfield_on[i].boundaryField()[4][j] - Tref));
		//stantonOLD(j,0) = stantonOLD(j,0)/ (Re * Pr);
		//stanton(j,0) = (gradient_t(j,0)*lambda*Pr) / (Ub *nu * (problem->Tfield_on[i].boundaryField()[4][j] - Tref));
		stanton(j,0) = (gradient_t(i,0)*nu) / (Ub *lambda*Pr* (problem->Tfield_on[i].boundaryField()[4][j] - Tref));
		CoordHeater(j,0) = problem->Tfield_on[i].mesh().boundaryMesh()[4].faceCentres()[j].x();
		Cf(j,0) = (2*nu*-BC5[j].component(0))/(Ub*Ub)*1000;

	   }

	   ITHACAstream::exportMatrix(stanton, "stantonFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   //ITHACAstream::exportMatrix(stantonOLD, "stantonFOMOLD", "eigen",
                  //             "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(CoordHeater, "CoordHeaterFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	   ITHACAstream::exportMatrix(Cf, "CfFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));


	   label stantonLWSize = problem->Tfield_on[i].boundaryFieldRef()[3].size();
	   Eigen::MatrixXd stantonLW(stantonLWSize, 1);
	   //Eigen::MatrixXd stantonLWOLD(stantonLWSize, 1);
	   Eigen::MatrixXd CoordLW(stantonSize, 1);
	   Eigen::MatrixXd CfLW(stantonSize, 1);
	   vectorField BC5LW = problem->Ufield_on[i].boundaryField()[3].snGrad();

	   for (label j = 0; j < stantonLWSize; j++)
	   {
	   	//stantonLWOLD(j,0) = (gradient_t(i,0) * h) / (lambda * (problem->Tfield_on[i].boundaryField()[3][j] - Tref));
		//stantonLWOLD(j,0) = stantonLWOLD(j,0)/ (Re * Pr);
	        //stantonLW(j,0) = (gradient_t(j,0)*lambda*Pr) / (Ub *nu * (problem->Tfield_on[i].boundaryField()[3][j] - Tref));
		stantonLW(j,0) = (gradient_t(i,0)*nu) / (Ub *lambda*Pr* (problem->Tfield_on[i].boundaryField()[3][j] - Tref));
		
		CoordLW(j,0) = problem->Tfield_on[i].mesh().boundaryMesh()[3].faceCentres()[j].x();
		CfLW(j,0) = (2*nu*-BC5LW[j].component(0))/(Ub*Ub)*1000;
	   }


	  ITHACAstream::exportMatrix(stantonLW, "stantonLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	  //ITHACAstream::exportMatrix(stantonLWOLD, "stantonLWFOMOLD", "eigen",
                          //     "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(CoordLW, "CoordLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

	  ITHACAstream::exportMatrix(CfLW, "CfLWFOM", "eigen",
                               "./ITHACAoutput/stanton/"+name(i+1));

    }

}


// ************************************************************************* //


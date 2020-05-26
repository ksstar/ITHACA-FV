	
	// As usual :  

	dimensionedScalar dt_fake
        (
            "dt_fake",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            scalar(1.0)
        );


	volVectorField U_rec("U_rec", Umodes[0]*0);
	for (label j = 0; j < Nphi_u; j++)
        {
            U_rec += Umodes[j] * a_o(j);
        }

	volVectorField U1aux("U1aux",problem->_U);

	U1aux =   dt*dt_fake * (-fvc::div(fvc::flux(U_rec),U_rec));
	volScalarField rhs = (1/dt)*(1/dt_fake)*fvc::div(U1aux);




	// InternalField (U_rec2) + Field containing the contribution of the boundary (U0) :

	volVectorField U_rec2("U_rec2", problem->_Ub); //Ub is zeroField with zeroGradient outlet and 0 inlet
	for (label j = 0; j < Nphi_u; j++)
        {
            U_rec2 += Umodes[j] * a_o(j);
        }

	
	volVectorField U_b("U_b", problem->_Ub);
	volVectorField U2aux("U2aux",problem->_U());

	U2aux =   dt*dt_fake * (
				-fvc::div(fvc::flux(U_rec2),U_rec2)
				-fvc::div(fvc::flux(U0),U0)
				-fvc::div(fvc::flux(U0),U_rec2)
				-fvc::div(fvc::flux(U_rec2),U0)

				+fvc::div(fvc::flux(U_b),U_rec2)
				+fvc::div(fvc::flux(U_rec2),U_b)
				+fvc::div(fvc::flux(U_b),U0)
				+fvc::div(fvc::flux(U0),U_b)
				);
	
	volScalarField rhs2 = (1/dt)*(1/dt_fake)*fvc::div(U2aux);


/////////////


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


	volVectorField U_rec("U_rec", problem->_Ub);
	for (label j = 0; j < Nphi_u; j++)
        {
            U_rec += Umodes[j] * a_o(j);
        }

	Vector<double> inl(1, 0, 0);
    	//ITHACAutilities::changeBCtype(U_rec,"fixedValue",0);
    	ITHACAutilities::assignBC(U_rec,0,inl);

	volVectorField Faux(problem->_U);

	Faux = dt * dt_fake*(nu * fvc::laplacian(nu_fake,U_rec)-fvc::div(fvc::flux(U_rec),U_rec));
	
	
	volScalarField p_b("p",problem->_p());
	    fvScalarMatrix pEqn = fvm::laplacian(p_b);
	    pEqn.setReference(0, 0);
	    solve(pEqn ==  (1/dt)*(1/dt_fake)*fvc::div(Faux)); 

	//volVectorField U_rec("U_rec", problem->Ufield[i]);
	///ITHACAstream::exportSolution(U_rec , name(i), "./ITHACAoutput/intersol");
	
	//volVectorField U1aux("U1aux",problem->_U);
	//U1aux =   dt*dt_fake * (-fvc::div(fvc::flux(U_rec),U_rec));
	//volScalarField rhs = (1/dt)*(1/dt_fake)*fvc::div(U1aux);

	/*ITHACAstream::exportSolution(U1aux , name(i), "./ITHACAoutput/intersol");
	//ITHACAstream::exportSolution(rhs , name(i), "./ITHACAoutput/intersol");

	volVectorField U_rec2("U_rec2", problem->_Ub);
	for (label j = 0; j < Nphi_u; j++)
        {
            U_rec2 += Umodes[j] * a_o(j);
        }
	//volVectorField U_rec2("U_rec", problem->Ufield[i]);
	Vector<double> v0(0, 0, 0);
	//ITHACAutilities::changeBCtype( U_rec2,"fixedValue",1);
	//ITHACAutilities::assignBC( U_rec2,1,v0);
	ITHACAstream::exportSolution(U_rec2 , name(i), "./ITHACAoutput/intersol");

	volVectorField U0("U0",problem->_Ub());
	ITHACAutilities::assignBC( U0,0,inl);
	ITHACAutilities::changeBCtype( U0,"fixedValue",1);
	ITHACAutilities::assignBC( U0,1,v0);
	ITHACAstream::exportSolution(U0 , name(i), "./ITHACAoutput/intersol");

	volVectorField U_b("U_b", problem->_Ub);

	volVectorField U2aux("U2aux",problem->_U());
	U2aux =   dt*dt_fake * (
				-fvc::div(fvc::flux(U_rec2),U_rec2)
				-fvc::div(fvc::flux(U0),U0)
				-fvc::div(fvc::flux(U0),U_rec2)
				-fvc::div(fvc::flux(U_rec2),U0)
				//+fvc::div(fvc::flux(U_b),U_rec2)
				//+fvc::div(fvc::flux(U_rec2),U_b)
				//+fvc::div(fvc::flux(U_b),U0)
				//+fvc::div(fvc::flux(U0),U_b)
				);

	ITHACAstream::exportSolution(U2aux , name(i), "./ITHACAoutput/intersol");

  	volVectorField Uout("Uout", problem->_U());
	Uout =  U2aux-U1aux;
	ITHACAstream::exportSolution(Uout , name(i), "./ITHACAoutput/intersol");

	volScalarField rhs2 = (1/dt)*(1/dt_fake)*fvc::div(U2aux);
	//ITHACAstream::exportSolution(rhs2 , name(i), "./ITHACAoutput/intersol");*/

	
	b_p = ITHACAutilities::get_coeffs(p_b,Pmodes);
	//std::cout << "  b_p   " <<  b_p << std::endl;
 	//b_q = ITHACAutilities::get_coeffs(rhs2,Pmodes);
	//std::cout << "  b_q  " <<  b_q << std::endl;

        List<Eigen::MatrixXd> RedLinSysP = problem->RedLinSysP;
        RedLinSysP[1] =(1/dt)*RedLinSysP[1] +
			(1/dt) * M4 	   + b_p;
 	
        b_n = reducedProblem::solveLinearSysAxb(RedLinSysP, xx, xx, presidual);


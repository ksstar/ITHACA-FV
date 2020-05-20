	
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


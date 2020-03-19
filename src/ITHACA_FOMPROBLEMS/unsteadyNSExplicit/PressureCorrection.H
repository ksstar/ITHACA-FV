U.correctBoundaryConditions();

solve(pEqn ==  (1/dt)*fvc::div(U));

U = U - dt*fvc::grad(p);

U.correctBoundaryConditions();

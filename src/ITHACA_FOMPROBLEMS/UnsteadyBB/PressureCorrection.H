 U.correctBoundaryConditions();

solve(pEqn == fvc::div(U));

U = U - fvc::grad(p);

U.correctBoundaryConditions();

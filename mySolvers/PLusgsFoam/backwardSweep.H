{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = 0;
        vector rRhoU = vector(0,0,0);
        scalar rRhoE = 0;
        scalar fRho  = 0;
        vector fRhoU = vector(0,0,0);
        scalar fRhoE = 0;

        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index
            if (mesh.isInternalFace(faceI))
            {
                label nei = neighbour[faceI];
                label own = owner[faceI];

                // we use the mesh ordering property: own < nei
                if (nei != cellI)  // then data[nei] belongs to upper triangle
                {
                    scalar dvol = mag( (mesh.C()[own] - mesh.C()[nei]) & Sf[faceI] );
                    scalar eps = min(1.0, max(magSqr(U[nei]/a[nei]), 0.04));
                    scalar ac = 0.25 * lusgsOmega * ( (1.0+eps)*mag((U[nei] & Sf[faceI]) - meshPhi[faceI])  
					+ Foam::sqrt(magSqr((eps-1.0)*(U[nei]& Sf[faceI])) + eps*magSqr(2.0*a[nei]*magSf[faceI])));

                    scalar av = sqr(magSf[faceI]) / dvol * nuMax[nei];
                    scalar rho1  = rho[nei] + dRho[nei];
                    vector rhoU1 = rhoU[nei] + dRhoU[nei];
                    scalar rhoE1 = rhoE[nei] + dRhoE[nei];
                    scalar p1    = (kappa[nei]-1) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];
                    scalar phi0r = phi0 - meshPhi[faceI];
                    scalar phi1r = phi1 - meshPhi[faceI];
                    
                    rRho  += (ac + av) * dRho[nei];
                    fRho  -= 0.5 * (rho1*phi1r - rho[nei]*phi0r); 
                    rRhoU += (ac + av) * dRhoU[nei]; 
                    fRhoU -= 0.5 * ( (rhoU1*phi1r - rhoU[nei]*phi0r) + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += (ac + av) * dRhoE[nei];
                    fRhoE -= 0.5 * ((rhoE1*phi1r+p1*phi1) - (rhoE[nei]*phi0r+p[nei]*phi0)); 
                }
            }
        }
        scalar eps = min(1.0, max(magSqr(U[cellI]/a[cellI]), 0.04));
        scalar dp =(1 - eps) / magSqr(a[cellI]) * (kappa[cellI]-1.0) *
                   (0.5*magSqr(U[cellI])*fRho- (U[cellI] & fRhoU) + fRhoE);

        dRho[cellI]  = dRhoStar[cellI]  + (rRho + fRho  - dp)/ D[cellI];
        dRhoU[cellI] = dRhoUStar[cellI] + (rRhoU+ fRhoU - dp*U[cellI]) / D[cellI];
        dRhoE[cellI] = dRhoEStar[cellI] + (rRhoE+ fRhoE - dp*(h[cellI]+0.5*magSqr(U[cellI])))/ D[cellI];
    }
}

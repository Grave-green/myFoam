{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
	    const scalarField& pMeshPhi = meshPhi.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];

            scalar eps = min(1.0, max(magSqr(U[ci]/a[ci]), 0.04));
             
            scalar ac = 0.25 * lusgsOmega * ( (1.0+eps)*mag((U[ci] & pSf[facei]) - pMeshPhi[facei])  
					+ Foam::sqrt(magSqr((eps-1.0)*U[ci]& pSf[facei]) + eps*magSqr(2.0*a[ci]*mag(pSf[facei]))));
            D[ci] += ac;
        }
    }

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar fRho  = mesh.V()[cellI]/dt.value() * dRho[cellI];
        vector fRhoU = mesh.V()[cellI]/dt.value() * dRhoU[cellI];
        scalar fRhoE = mesh.V()[cellI]/dt.value() * dRhoE[cellI];

        scalar rRho  = 0.0;
        vector rRhoU = vector(0,0,0);
        scalar rRhoE = 0.0;

        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index
            
            if (mesh.isInternalFace(faceI))
            {
                // Calculate D
                label own = owner[faceI];
                label nei = neighbour[faceI];
	        	scalar dvol = mag( (mesh.C()[own] - mesh.C()[nei]) & Sf[faceI] );

                scalar eps = min(1.0, max(magSqr(U[cellI]/a[cellI]), 0.04));
                scalar ac = 0.25 * lusgsOmega * ( (1.0+eps)*mag((U[cellI] & Sf[faceI]) - meshPhi[faceI])  
					+ Foam::sqrt(magSqr((eps-1.0)*(U[cellI]& Sf[faceI])) + eps*magSqr(2.0*a[cellI]*(magSf[faceI]))));
 
                scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];

                D[cellI] += ac + av;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    scalar eps = min(1.0, max(magSqr(0.5*U[own]/a[own]), 0.04));
                    scalar ac = 0.25 * lusgsOmega * ( (1.0+eps)*mag((U[own] & Sf[faceI]) - meshPhi[faceI])  
					+ Foam::sqrt(magSqr((eps-1.0)*(U[own]& Sf[faceI])) + eps*magSqr(2.0*a[own]*magSf[faceI])));
                    scalar av = sqr(magSf[faceI]) / dvol * nuMax[own];
                    
                    scalar rho1  = rho[own] + dRhoStar[own];
                    vector rhoU1 = rhoU[own] + dRhoUStar[own];
                    scalar rhoE1 = rhoE[own] + dRhoEStar[own];
                    scalar p1    = (kappa[own]-1) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];
		            scalar phi0r = phi0 - meshPhi[faceI];
		            scalar phi1r = phi1 - meshPhi[faceI];

                    rRho  += (ac + av) * dRhoStar[own];
                    fRho  += 0.5 * (rho1*phi1r - rho[own]*phi0r); 
                    rRhoU += (ac + av) * dRhoUStar[own];
                    fRhoU += 0.5 * ( (rhoU1*phi1r - rhoU[own]*phi0r) + (p1 - p[own])*Sf[faceI] );
                    rRhoE += (ac + av) * dRhoEStar[own];
                    fRhoE += 0.5 * (rhoE1*phi1r+p1*phi1 - (rhoE[own]*phi0r+p[own]*phi0)); 
 
                }
            }
        }

        scalar eps = min(1.0, max(magSqr(U[cellI]/a[cellI]), 0.04));
        scalar dp =(1 -eps) / magSqr(a[cellI]) * (kappa[cellI]-1.0) *
                   (0.5*magSqr(U[cellI])*fRho- (U[cellI] & fRhoU) + fRhoE);

        dRhoStar[cellI]  = (mesh.V()[cellI]/dt.value()* ddtRho[cellI] + rRho  + fRho  - dp)/ D[cellI];
        dRhoUStar[cellI] = (mesh.V()[cellI]/dt.value()*ddtRhoU[cellI] + rRhoU + fRhoU - dp*U[cellI] ) / D[cellI];
        dRhoEStar[cellI] = (mesh.V()[cellI]/dt.value()*ddtRhoE[cellI] + rRhoE + fRhoE - dp*(h[cellI]+0.5*magSqr(U[cellI])) )/ D[cellI];
    }
}

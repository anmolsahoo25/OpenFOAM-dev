/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mprkFoam

Description
    Solves the advection problem using MPRK method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "dynamicRefineFvMesh.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "\nCreating the mesh as dynamicRefineFvMesh" << endl;

    Foam::dynamicRefineFvMesh mesh
    (
	    Foam::IOobject
	    (
		    Foam::fvMesh::defaultRegion,
		    runTime.timeName(),
		    runTime,
		    Foam::IOobject::MUST_READ
	    )
    );

    simpleControl simple(mesh);

    #include "createFields.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Marking the coarse and fine zones

        labelList fineCellRegion(mesh.nCells(), 0);
        labelList coarseCellRegion(mesh.nCells(), 0);
        
        scalar coarseVolume = 0.0f;
        scalar fineVolume = mesh.V()[0];

        // Finding the coarse volume and fine volumes values
        forAll(mesh.C(), i) {
            if(mesh.V()[i] > coarseVolume) {
                coarseVolume = mesh.V()[i];
            }
            if(mesh.V()[i] < fineVolume) {
                fineVolume = mesh.V()[i];
            }
        }

        // Marking the regions
        forAll(mesh.C(), i) {
            if(cmptMag(mesh.V()[i] - coarseVolume) > 1e-10) {
                fineCellRegion[i] = 2;
            }
            else {
                coarseCellRegion[i] = 2;
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Creating fast-slow interface for stability

        // Number of cells around cell zone to create as interface
        int numLayers = 5;

        for(int i = 0; i < numLayers; i++) {
            forAll(mesh.C(), j) {
                if(fineCellRegion[j] != 2) {
                    forAll(mesh.cellCells()[j], k) {
                        if(fineCellRegion[mesh.cellCells()[j][k]] == 2) {
                            fineCellRegion[j] = 2;
                        }
                    }
                }

                if(coarseCellRegion[j] != 2) {
                    forAll(mesh.cellCells()[j], k) {
                        if(coarseCellRegion[mesh.cellCells()[j][k]] == 2) {
                            coarseCellRegion[j] = 2;
                        }
                    }
                }
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Handling cyclic patches
        // Make sure that the any cells adjacent to cyclic faces are in the subset
        
        const fvBoundaryMesh& boundary = mesh.boundary();
        forAll(boundary, patchI) {
            if(boundary[patchI].coupled() == true) {
                forAll(boundary[patchI].faceCells(), cellI) {
                    label cellIndex = boundary[patchI].faceCells()[cellI];
                    fineCellRegion[cellIndex] = 2;
                    coarseCellRegion[cellIndex] = 2;
                    forAll(mesh.cellCells()[cellIndex], k) {
                        label n1 = mesh.cellCells()[cellIndex][k];
                        fineCellRegion[n1] = 2;
                        coarseCellRegion[n1] = 2;
                        forAll(mesh.cellCells()[n1], l) {
                            label n2 = mesh.cellCells()[n1][l];
                            fineCellRegion[n2] = 2;
                            coarseCellRegion[n2] = 2;
                        }
                    }
                }
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Creating the mesh subsets
        
        fvMeshSubset subsetFine(mesh);
        fvMeshSubset subsetCoarse(mesh);

        subsetFine.setLargeCellSubset(fineCellRegion, 2);
        subsetCoarse.setLargeCellSubset(coarseCellRegion, 2);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Extracting the fields for the different zones

        // Creating the fields for the fine zone
        volScalarField cf = subsetFine.interpolate(c)();
        surfaceScalarField phif = subsetFine.interpolate(phi)();

        // Creating the fields for the coarse zone
        volScalarField cc = subsetCoarse.interpolate(c)();
        surfaceScalarField phic = subsetCoarse.interpolate(phi)();
        
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Performing one time step increment for the various zones

        scalar scaleFactor = fineVolume / coarseVolume;
        int subSteps = static_cast<int>(1/scaleFactor);

        dimensionedScalar deltaT1(
            "deltaT1",
            dimensionSet(0,0,1,0,0,0,0),
            1.0
        );

        dimensionedScalar deltaTf(
            "deltaTf",
            dimensionSet(0,0,1,0,0,0,0),
            scaleFactor*runTime.deltaTValue()
        );

        dimensionedScalar deltaTf2(
            "deltaTf2",
            dimensionSet(0,0,1,0,0,0,0),
            0.5*scaleFactor*runTime.deltaTValue()
        );

        dimensionedScalar deltaT(
            "deltaT",
            dimensionSet(0,0,1,0,0,0,0),
            runTime.deltaTValue()
        );

        // Time stepping for the fine zone
        volScalarField cf_temp(cf);

        volScalarField RHSf(cf);
        forAll(RHSf, i) {
            RHSf[i] = 0.0;
        }

        // Substepping for the fine region
        for(int i = 0; i < subSteps; i++) {
            volScalarField cf_temp1("cf_temp1", cf_temp);
            volScalarField k1f = fvc::div(phif, cf_temp);
            cf_temp = cf_temp + k1f*deltaTf;
            volScalarField k2f = fvc::div(phif, cf_temp);
            RHSf = RHSf - 0.5*scaleFactor*runTime.deltaTValue()*deltaT1*(k1f + k2f);
            cf_temp = cf_temp + deltaTf2*(k1f + k2f);
        }

        // Time stepping for the coarse zone
        volScalarField cc_temp(
            IOobject(
                "cc_temp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            cc
        );

        volScalarField RHSc(cc);
        forAll(RHSc, i) {
            RHSc[i] = 0.0;
        }

        // Substepping for the fine region
        volScalarField k1c = fvc::div(phic, cc_temp);
        cc_temp = cc_temp + k1c*deltaT;
        volScalarField k2c = fvc::div(phif, cc_temp);
        RHSc = RHSc - 0.5*runTime.deltaTValue()*deltaT1*(k1c + k2c);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Remapping the RHSf and RHSc obtained to the main RHS
        
        volScalarField RHS(c);
        forAll(RHS, i) {
            RHS[i] = 0.0;
        }

        // If the cell is coarse, take value from coarse subset
        forAll(subsetCoarse.cellMap(), i) {
            label cellIndex = subsetFine.cellMap()[i];
            if(cmptMag(mesh.V()[cellIndex] - coarseVolume) < 1e-10) {
                RHS[cellIndex] = RHSc[i];
            }
        }

        // If the cell is fine, take value from fine subset and also set neighbours to fine for stability zone
        forAll(subsetFine.cellMap(), i) {
            label cellIndex = subsetFine.cellMap()[i];
            if(cmptMag(mesh.V()[cellIndex] - coarseVolume) > 1e-10) {
                RHS[cellIndex] = RHSf[i];
            }
            else {
                forAll(mesh.cellCells()[cellIndex], k) {
                    if(cmptMag(mesh.V()[mesh.cellCells()[cellIndex][k]] - coarseVolume) > 1e-10) {
                        RHS[cellIndex] = RHSf[i];
                    }
                }
            }
        }


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Solve the equation for c
        
        fvScalarMatrix cEqn(fvm::ddt(c) == RHS);
        cEqn.solve();

        runTime.write();

        mesh.update();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

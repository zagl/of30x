/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    thinWalls

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "hexRef8.H"
#include "cellSet.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "topoSetSource.H"
#include "IOobjectList.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createNamedPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");

    const word dictName("snappyHexMeshDict");
#   include "setSystemMeshDictionaryIO.H"
    IOdictionary dict(dictIO);

    const dictionary& geometryDict = dict.subDict("geometry");
    const dictionary& surfacesDict = 
        dict.subDict("castellatedMeshControls").subDict("refinementSurfaces");

    hexRef8 meshCutter(mesh);
    const vectorField& cellCentres = mesh.cellCentres();
    const cellList& cells = mesh.cells();
    const scalar baseLength = meshCutter.level0EdgeLength();
    const labelList& cellLevel = meshCutter.cellLevel(); 

    cellZoneMesh& cellZones = mesh.cellZones();
    
    PackedBoolList isAssigned(mesh.nCells(), false);

    labelHashSet facesInZones(mesh.nFaces()/10);

    forAllConstIter(dictionary, geometryDict, iter)
    {
        const fileName surfFileName = iter().keyword();
        const word surfName = iter().dict().lookup("name");
        const dictionary& surfaceDict = surfacesDict.subDict(surfName);

        if ( surfaceDict.found("cellZone") )
        {
            const triSurface surf(runTime.constantPath()/"triSurface"/surfFileName);
            const vectorField& normals = surf.faceNormals();
            const triSurfaceSearch querySurf(surf);
            const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

            const word faceZoneName = surfaceDict.lookup("faceZone");
            const word cellZoneName = surfaceDict.lookup("cellZone");

            cellZoneSet regionCellZone
            (
                mesh,
                cellZoneName,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            labelList& regionCells = regionCellZone.addressing();
            
            regionCells.clear();

            regionCellZone.updateSet();
            regionCellZone.write();

            faceZoneSet regionFaceZone
            (
                mesh,
                faceZoneName,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            forAll( cellCentres, cellI )
            {
                if ( !isAssigned[cellI] )
                {
                    const point& cellCentre = cellCentres[cellI];
                    
                    if ( tree.bb().contains(cellCentre) )
                    {
                        if (tree.getVolumeType(cellCentre) == volumeType::INSIDE)
                        {
                            regionCells.append(cellI);
                            isAssigned[cellI] = true;
                        }
                    }
                    else
                    {
                        scalar localLength = baseLength / pow(2, cellLevel[cellI]);
                        pointIndexHit pHit1 =
                            tree.findNearest(cellCentre, pow(localLength*1., 2));
                        if ( pHit1.hit() )
                        {
                            label triangleI = pHit1.index();
                            point nearestSurfPoint = pHit1.hitPoint();
                            vector localNormal = normals[triangleI];
                            point searchStart = nearestSurfPoint - localNormal*1e-6;
                            point searchEnd = searchStart - localNormal*(100*baseLength);
                            pointIndexHit pHit2 = tree.findLine(searchStart, searchEnd);
                            if ( pHit2.hit() )
                            {
                                point nextSurfPoint = pHit2.hitPoint();
                                scalar wallThickness = mag(nextSurfPoint - nearestSurfPoint);
                                if ( wallThickness <= localLength*1.7 )
                                {
                                    regionCells.append(cellI);
                                    isAssigned[cellI] = true;
                                }
                            }
                        }
                    }
                }
            }

            regionCellZone.updateSet();
            regionCellZone.write();

            labelList nFaceCells(mesh.nFaces(), 0);
            forAll( regionCells, i )
            {
                label cellI = regionCells[i];
                const cell& faces = cells[cellI];
                forAll( faces, i )
                {
                    nFaceCells[faces[i]]++;
                }
            }

            DynamicList<label> regionFaces;
            forAll( nFaceCells, faceI )
            {
                if ( nFaceCells[faceI] == 1 && ! facesInZones[faceI] )
                {
                    regionFaces.append(faceI);
                    facesInZones.insert(faceI);
                }
            }
            regionFaceZone.addressing().transfer(regionFaces);
            regionFaceZone.updateSet();
            regionFaceZone.write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

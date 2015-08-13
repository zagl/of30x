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
    newMoveCells

Description
    Utility to move cell centres into thin triSurface walls. Meant as 
    preprocessing step for snappyHexMesh.

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
#include "plane.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    argList::validArgs.append("input surfaceFile");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");
    const fileName surfName = args[1];

    triSurface surf(runTime.constantPath()/"triSurface"/surfName);
    const vectorField& normals = surf.faceNormals();

    const pointField& points = mesh.points();
    const labelListList& cellPoints = mesh.cellPoints();
    labelList edgeLabels(mesh.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    const edgeList edges = mesh.edges();

    meshSearch queryMesh(mesh);
    triSurfaceSearch querySurf(surf);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

    Info << "Find points near triSurface" << endl << endl;

    pointField newLocations(mesh.nPoints());
    DynamicList<label> movedPoints(mesh.nPoints()/4);
    PackedBoolList pointsNotMoved(mesh.nPoints(), true);
    PackedBoolList cellsNotMoved(mesh.nCells(), true);
    PackedBoolList cellsPartlyMoved(mesh.nCells(), false);

    //DynamicList<point> movePoints(mesh.nPoints());

    hexRef8 meshCutter(mesh);
    const vectorField& cellCentres = mesh.cellCentres();
    const labelList& cellLevel = meshCutter.cellLevel(); 
    const scalar baseLength = meshCutter.level0EdgeLength();

    label outside = 0;
    forAll(cellCentres, i)
    {
        vector centre = cellCentres[i];
        if ( tree.getVolumeType(centre) == volumeType::OUTSIDE )
        {
            scalar localLength = baseLength / pow(2, cellLevel[i]);
            pointIndexHit pHit =
                tree.findNearest(centre, pow(localLength*0.5, 2));
            if ( pHit.hit() )
            {
                label triangleI = pHit.index();
                point p0 = pHit.hitPoint();
                vector normal = normals[triangleI];
                point pStart = p0 - normal*1e-6;
                point pEnd = pStart - normal*(100*baseLength);
                pointIndexHit pNext = tree.findLine(pStart, pEnd);
                if ( pNext.hit() )
                {
                    point p1 = pNext.hitPoint();
                    scalar distance = mag(p0 - p1);
                    if ( distance <= localLength )
                    {
                        point pCentre = (p0 + p1)/2;
                        label changeCell = queryMesh.findNearestCell(pCentre, i);
                        point changeCentre = cellCentres[changeCell];

                        if (
                            cellsNotMoved[changeCell]
                            && tree.getVolumeType(changeCentre)
                               == volumeType::OUTSIDE
                        )
                        {
                            scalar changeLength = 
                                baseLength / pow(2, cellLevel[changeCell]);
                            labelList cPoints = cellPoints[changeCell];
                            pointIndexHit pHit =
                                tree.findNearest(changeCentre, pow(changeLength, 2));
                            point changeNearest = pHit.hitPoint();
                            scalar changeDist = mag(changeNearest-changeCentre);
                            vector changeNormal =
                                (changeNearest - changeCentre)/changeDist;
                            point newCentre =
                                changeNearest + changeNormal*distance*0.2;

                            vector moveVector = newCentre - changeCentre;
                            plane triPlane(newCentre, changeNormal);

                            forAll( cPoints, i )
                            {
                                label pointI = cPoints[i];
                                label oppPointI = cPoints[i^6];
                                if ( pointsNotMoved[pointI] )
                                {
                                    point final;
                                    //if ( pointsNotMoved[oppPointI] )
                                    //{
                                        point moved = points[pointI] + moveVector;
                                        point near = triPlane.nearestPoint(moved);
                                        final = moved + (near - moved)*0.2;
                                    //}
                                    //else
                                    //{
                                    //    point oppPoint = newLocations[oppPointI];
                                    //    final = 2*newCentre - oppPoint;
                                    //    Info<< changeCell << " ";
                                    //}
                                    movedPoints.append(pointI);
                                    newLocations[pointI] = final;
                                    pointsNotMoved[pointI] = false;
                                }
                                else
                                {
                                    cellsPartlyMoved[changeCell] = true;
                                }
                            }
                            cellsNotMoved[changeCell] = false;
                        }
                    }
                }
            }
        }
    }

    polyTopoChange meshMod(mesh);
    forAll(movedPoints, i)
    {
        label pointI = movedPoints[i];
        point newLocation = newLocations[pointI];

        meshMod.modifyPoint(pointI, newLocation, -1, true);
    }
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    const vectorField& newCellCentres = mesh.cellCentres();

    forAll( cellsPartlyMoved, cellI )
    {
        if ( cellsPartlyMoved[cellI] )
        {
            if ( tree.getVolumeType(newCellCentres[cellI]) == volumeType::OUTSIDE )
                Info<< " " << cellI;
        }
    }


    Info<< nl << nl;

    forAll( cellsPartlyMoved, cellI )
    {
        if ( cellsPartlyMoved[cellI] )
        {
            if ( tree.getVolumeType(newCellCentres[cellI]) == volumeType::INSIDE )
                Info<< " " << cellI;
        }
    }




    if (!overwrite)
    {
        runTime++;
    }
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing morphMesh to time " << runTime.timeName() << endl << endl;

    mesh.write(); 

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

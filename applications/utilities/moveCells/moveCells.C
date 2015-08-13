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
    moveMeshPoints

Description
    Utility to move Points of a Mesh onto a near triSurface.

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
//    argList::noParallel();
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

    DynamicList<point> newLocations(mesh.nPoints());

    List<List<point> > movePoints(mesh.nPoints());

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
                        labelList cPoints = cellPoints[changeCell];
                        point changeCentre = cellCentres[changeCell];

                        vector moveVector = pCentre - changeCentre;
                        //vector moveVector = pStart - centre;
                        plane triPlane(pCentre, normal);
                        forAll( cPoints, i )
                        {
                            label pointI = cPoints[i];
                            point moved = points[pointI] + moveVector;
                            point near = triPlane.nearestPoint(moved);
                            point final = moved + (near - moved)*0.5;
                            vector vFinal = final - points[pointI];

                            //movePoints[pointI].append(vFinal);
                            movePoints[pointI].append(moveVector);
                        }
                    }
                }
            }
        }
    }


    polyTopoChange meshMod(mesh);
    forAll(movePoints, i)
    {
        List<point>& newPoints = movePoints[i];
        label number = newPoints.size();
        if ( number > 0 )
        {
            scalar minMag = GREAT;
            label smallest;
            forAll( newPoints, pointI )
            {
                vector v = newPoints[pointI];
                if( mag(v) < minMag )
                {
                    smallest = pointI;
                }
            }
            //point newLocation = points[i] + newPoints[smallest];
            point newLocation = vector::zero;
            forAll( newPoints, pointI )
            {
                point newPoint = newPoints[pointI] + points[i];
                newLocation += newPoint;
            }
            newLocation /= number;
            meshMod.modifyPoint(i, newLocation, -1, true);
        }
    }
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);



//    forAll(edgeLabels, i)
//    {
//
//
//        pointIndexHit pHit = tree.findNearest(pStart, eTol);
//
//
//
//
//        label edgeI = edgeLabels[i];
//        const edge e = edges[edgeI];
//        const point pStart = points[e.start()] ;
//        const point pEnd = points[e.end()] ;
//        const vector eVec(pEnd - pStart);
//        const scalar eMag = mag(eVec);
//        const vector n(eVec/(eMag + VSMALL));
//        const point tolVec = 1e-6*eVec;
//        const scalar eTol = tol * eMag;
//        point p0 = pStart - tolVec;
//        const point p1 = pEnd + tolVec;
//        bool foundStart = false;
//        bool foundEnd = false;
//        pointIndexHit nearHitStart;
//        pointIndexHit nearHitEnd;
//        while(true)
//        {
//            pointIndexHit pHit = tree.findLine(p0, p1);
//
//            if(pHit.hit())
//            {
//                if (mag(pHit.hitPoint() - pStart) < eTol && !foundStart)
//                {
//                    foundStart = true;
//                }
//                else if (mag(pHit.hitPoint() - pEnd) < eTol)
//                {
//                    foundEnd = true;
//
//                    if (mag(pHit.hitPoint() - pEnd) < 1e-6*eMag)
//                    {
//                        // Reached end.
//                        break;
//                    }
//                }
//                p0 = pHit.hitPoint() + tolVec;
//            }
//            else
//            {
//                // No hit.
//                break;
//            }
//        }
//
//        if(foundStart)
//        {
//            pointIndexHit pHit = tree.findNearest(pStart, eTol);
//            if (pHit.hit())
//            {
//                movePoints.append(e.start());
//                newLocations.append(pHit.hitPoint());
//            }
//        }
//        if(foundEnd)
//        {
//            pointIndexHit pHit = tree.findNearest(pEnd, eTol);
//            if (pHit.hit())
//            {
//                movePoints.append(e.end());
//                newLocations.append(pHit.hitPoint());
//            }
//        }
//    }
//
//    movePoints.shrink();
//    newLocations.shrink();

    //Info<< "Moving " << movePoints.size() << " points"  << endl << endl;

    //polyTopoChange meshMod(mesh);
    //forAll(movePoints, i)
    //{
    //    label movePointI = movePoints[i];
    //    point newLocationI = newLocations[i];

    //    meshMod.modifyPoint(movePointI, newLocationI, -1, true);
    //}
    //autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

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

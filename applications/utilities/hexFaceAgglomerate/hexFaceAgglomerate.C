/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    hexFaceAgglomerate

Description
    Agglomerate boundary faces by hex regions.
    It writes a map from the fine to coarse grid.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "boundaryMesh.H"
#include "boundBox.H"
#include "volFields.H"
#include "CompactListList.H"
#include "unitConversion.H"
#include "labelListIOList.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "labelVector.H"
#include "PatchTools.H"
#include "uindirectPrimitivePatch.H"
#include "OFstream.H"
#include "meshTools.H"
#include "ListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("viewFactorsDict");

    #include "setConstantMeshDictionaryIO.H"

    // Read control dictionary
    const IOdictionary agglomDict(dictIO);

    bool writeAgglom = readBool(agglomDict.lookup("writeFacesAgglomeration"));


    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    OFstream debug("debug.obj");


    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        boundary.size()
    );

    label nCoarseFaces = 0;

    forAllConstIter(dictionary, agglomDict, iter)
    {
        labelList patchIds = boundary.findIndices(iter().keyword());


        forAll(patchIds, i)
        {
            label patchI =  patchIds[i];
            const polyPatch& pp = boundary[patchI];

            if (!pp.coupled())
            {
                Info << "\nAgglomerating patch : " << pp.name() << endl;

                dictionary controlDict = agglomDict.subDict(pp.name());
                label nFacesInCoarsestLevel = readLabel(controlDict.lookup("nFacesInCoarsestLevel"));
                scalar featureAngle = controlDict.lookupOrDefault<scalar>("featureAngle", 0);
                const scalar minCos = Foam::cos(degToRad(featureAngle));

                boundBox bb(pp.localPoints());

                label maxOneDir = nFacesInCoarsestLevel;
                scalar stepLength = cmptMax( bb.span() ) / maxOneDir;

                point bbMin = bb.min() - vector(stepLength*.5, stepLength*.5, stepLength*.5);

                label nX = int( bb.span().x() / stepLength ) + 2;
                label nY = int( bb.span().y() / stepLength ) + 2;
                label nZ = int( bb.span().z() / stepLength ) + 2;


                labelListList agglomCellsFaces(nX*nY*nZ);


                const pointField& faceCentres = pp.faceCentres();
                const vectorField& faceNormals = pp.faceNormals();
                const labelListList& faceEdges = pp.faceEdges();
                const labelListList& faceFaces = pp.faceFaces();
                const labelListList& edgeFaces = pp.edgeFaces();
                const labelListList& pointEdges = pp.pointEdges();
                const edgeList& edges = pp.edges();
                const pointField& points = pp.points();

                labelList patchAgglomeration(pp.size());

                forAll(faceCentres, faceI)
                {
                    const point& faceCentre = faceCentres[faceI];

                    vector faceOffset = faceCentre - bbMin;

                    label xI = int( faceOffset.x() / stepLength );
                    label yI = int( faceOffset.y() / stepLength );
                    label zI = int( faceOffset.z() / stepLength );

                    label agglomCellI = xI + yI*nX + zI*nX*nY;

                    agglomCellsFaces[agglomCellI].append(faceI);

                }

                label coarseFaceI = 0;

                forAll( agglomCellsFaces, agglomCellI )
                {
                    labelList agglomCellFaces = agglomCellsFaces[agglomCellI];

                    if (agglomCellFaces.size() > 0)
                    {

                        labelListList edgesFaces(pp.nEdges());

                        forAll( agglomCellFaces, i )
                        {
                            label faceI = agglomCellFaces[i];
                            labelList faceIEdges = faceEdges[faceI];
                            forAll( faceIEdges, j )
                            {
                                label edgeI = faceIEdges[j];
                                edgesFaces[edgeI].append(faceI);
                            }
                        }

                        boolList borderEdge(pp.nEdges(), false);
                        DynamicList<label> borderEdges(pp.nEdges());

                        forAll( edgesFaces, edgeI )
                        {
                            labelList edgeFaces = edgesFaces[edgeI];

                            bool isBorderEdge = false;

                            if (edgeFaces.size() == 2)
                            {
                                if ((faceNormals[edgeFaces[0]] & faceNormals[edgeFaces[1]]) < minCos)
                                {
                                    isBorderEdge = true;
                                }
                            }
                            else
                            {
                                isBorderEdge = true;
                            }

                            if (isBorderEdge)
                            {
                                borderEdge[edgeI] = true;
                                borderEdges.append(edgeI);
                            }
                        }

                        label currentZone = 0;
                        labelList faceZone(pp.size(), -1);

                        forAll( agglomCellFaces, i )
                        {
                            label faceI = agglomCellFaces[i];

                            if (faceZone[faceI] == -1)
                            {
                                PatchTools::markZone(
                                    pp,
                                    borderEdge,
                                    faceI,
                                    currentZone,
                                    faceZone
                                );
                                faceZone[faceI] = currentZone;
                                currentZone++;
                            }

                            patchAgglomeration[faceI] = coarseFaceI + faceZone[faceI];
                        }

                        labelListList zones(currentZone);

                        forAll(faceZone, faceI)
                        {
                            if (faceZone[faceI] != -1)
                            {
                                zones[faceZone[faceI]].append(faceI);
                            }
                        }

                        forAll(zones, zoneI)
                        {
                            const labelList& zoneFaces = zones[zoneI];
                            uindirectPrimitivePatch upp
                            (
                                UIndirectList<face>(pp, zoneFaces),
                                pp.points()
                            );

                            if (upp.edgeLoops().size() != 1)
                            {
                                label outerLoop = -1;
                                label maxLoopPoints = -1;

                                forAll( upp.edgeLoops(), loopI)
                                {
                                    const labelList& edgeLoop = upp.edgeLoops()[loopI];
                                    if (edgeLoop.size() > maxLoopPoints)
                                    {
                                        outerLoop = loopI;
                                        maxLoopPoints = edgeLoop.size();
                                    }
                                }

                                List<List<scalar> > splits(3);
                                List<label> nSplits(3, 1);
                                forAll( upp.edgeLoops(), loopI)
                                {
                                    if (loopI != outerLoop)
                                    {
                                        const labelList& edgeLoop = upp.edgeLoops()[loopI];
                                        List<point> loopPoints;

                                        forAll( edgeLoop, i)
                                        {
                                            loopPoints.append(upp.localPoints()[edgeLoop[i]]);
                                            meshTools::writeOBJ(debug, upp.localPoints()[edgeLoop[i]]);
                                        }

                                        boundBox loopbb(loopPoints);
                                        label longestEdge = findMax(loopbb.span());
                                        splits[longestEdge].append(loopbb.midpoint()[longestEdge]);
                                        nSplits[longestEdge]++;
                                    }
                                }

                                labelListList order(3);
                                forAll(splits, i)
                                {
                                    splits[i].append(VGREAT);
                                    sortedOrder(splits[i], order[i]);
                                }

                                labelList subFaceZone(upp.size(), -1);
                                labelListList subFaceZones(nSplits[0]*nSplits[1]*nSplits[2]);
                                forAll( zoneFaces, i )
                                {
                                    const label& faceI = zoneFaces[i];
                                    const point& faceCentre = faceCentres[faceI];

                                    forAll(splits[2], zI)
                                    {
                                        scalar z = splits[2][order[2][zI]];
                                        forAll(splits[1], yI)
                                        {
                                            scalar y = splits[1][order[1][yI]];
                                            forAll(splits[0], xI)
                                            {
                                                scalar x = splits[0][order[0][xI]];

                                                if (subFaceZone[i] == -1)
                                                {
                                                    if (
                                                        faceCentre[0] < x
                                                     && faceCentre[1] < y
                                                     && faceCentre[2] < z
                                                    )
                                                    {
                                                        label agglomSubCellI = xI + yI*nSplits[0] + zI*nSplits[0]*nSplits[1];
                                                        subFaceZone[i] = agglomSubCellI;
                                                        subFaceZones[agglomSubCellI].append(faceI);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                forAll( subFaceZones, zoneI )
                                {
                                    const labelList& zoneFaces = subFaceZones[zoneI];
                                    if (zoneFaces.size() > 0)
                                    {
                                        forAll(zoneFaces, i)
                                        {
                                            label faceI = zoneFaces[i];
                                            faceZone[faceI] = currentZone;
                                            patchAgglomeration[faceI] = coarseFaceI + currentZone;

                                        }
                                        currentZone++;
                                    }
                                }









                            }

                        }


                        forAll( agglomCellFaces, i )
                        {
                            label faceI = agglomCellFaces[i];
                            patchAgglomeration[faceI] = coarseFaceI + faceZone[faceI];
                        }


                        coarseFaceI += currentZone;
                    }
                }



                finalAgglom[patchI] = patchAgglomeration;


                if (finalAgglom[patchI].size())
                {
                    nCoarseFaces += max(finalAgglom[patchI] + 1);
                }
            }
        }
    }


    // - All patches which are not agglomarated are identity for finalAgglom
    forAll(boundary, patchId)
    {
        if (finalAgglom[patchId].size() == 0)
        {
            finalAgglom[patchId] = identity(boundary[patchId].size());
        }
    }

    // Sync agglomeration across coupled patches
    labelList nbrAgglom(mesh.nFaces() - mesh.nInternalFaces(), -1);

    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled())
        {
            finalAgglom[patchId] = identity(pp.size());
            forAll(pp, i)
            {
                nbrAgglom[pp.start() - mesh.nInternalFaces() + i] =
                    finalAgglom[patchId][i];
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, nbrAgglom);
    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled() && !refCast<const coupledPolyPatch>(pp).owner())
        {
            forAll(pp, i)
            {
                finalAgglom[patchId][i] =
                    nbrAgglom[pp.start() - mesh.nInternalFaces() + i];
            }
        }
    }

    finalAgglom.write();

    if (writeAgglom)
    {
        globalIndex index(nCoarseFaces);
        volScalarField facesAgglomeration
        (
            IOobject
            (
                "facesAgglomeration",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("facesAgglomeration", dimless, 0)
        );

        label coarsePatchIndex = 0;
        forAll(boundary, patchId)
        {
            const polyPatch& pp = boundary[patchId];
            if (pp.size() > 0)
            {
                fvPatchScalarField& bFacesAgglomeration =
                    facesAgglomeration.boundaryField()[patchId];

                forAll(bFacesAgglomeration, j)
                {
                    bFacesAgglomeration[j] =
                        index.toGlobal
                        (
                            Pstream::myProcNo(),
                            finalAgglom[patchId][j] + coarsePatchIndex
                        );
                }

                coarsePatchIndex += max(finalAgglom[patchId]) + 1;
            }
        }

        Info<< "\nWriting facesAgglomeration" << endl;
        facesAgglomeration.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

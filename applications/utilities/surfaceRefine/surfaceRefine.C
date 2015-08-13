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
    surfaceRefine

Description
    Refine faces with edges longer than specified length.

\*---------------------------------------------------------------------------*/


#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");
    argList::validArgs.append("maximum length");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName = args[2];  
    const scalar maxLength = args.argRead<scalar>(3);

    Info<< "Reading surface from " << surfFileName << " ..." << endl << endl;

    const triSurface& surf1(surfFileName);
    triSurface surf2;
    triSurface surf = surf1;
    
    int loopCounter = 1;
    
    while(true)
    {
        const pointField points(surf.points());
        const labelList meshPoints(surf.meshPoints());
        const edgeList& edges = surf.edges();
        labelListList edgeFaces = surf.sortedEdgeFaces();
        
        List<bool> refineFaces(surf.size(), false);
        
        forAll(edges, i)
        {   
            const edge e = surf.edges()[i];
            const point& pStart = points[meshPoints[e.start()]];
            const point& pEnd = points[meshPoints[e.end()]];
            const vector eVec(pEnd - pStart);
            const scalar eMag = mag(eVec);
            
            if(eMag > maxLength)
            {
                labelList edgeFacesI = edgeFaces[i];
                
                forAll(edgeFacesI, i)
                {
                    label faceI = edgeFacesI[i];
                    refineFaces[faceI] = true;
                }
            }
        }
        
        DynamicList<label> refineF(surf.size());
        
        forAll(refineFaces, i)
        {
            if(refineFaces[i])
            {
                refineF.append(i);
            }
        } 

        if(refineF.size() > 0)
        {
            surf2 = triSurfaceTools::redGreenRefine(surf, refineF);
            surf = surf2;
        }
        else
        {
            surf2 = surf;
            break;
        }
        
        Info << "Loop " << loopCounter << ": " << nl
             << " Refined triangles: " << refineF.size() << nl;
             
        loopCounter++;
    }
    
    Info<< nl
        << "Original surface:" << endl
        << " triangles :" << surf1.size() << endl
        << " vertices(used):" << surf1.nPoints() << endl << endl
        << "Refined surface:" << endl
        << " triangles :" << surf2.size() << endl
        << " vertices(used):" << surf2.nPoints() << endl << endl;

    Info<< "Writing refined surface to " << outFileName << " ..." << endl;
    
    surf2.write(outFileName);
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

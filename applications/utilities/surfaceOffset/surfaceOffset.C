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
    surfaceOffset

Description
    Move points of surface away from original surface in normals direction.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Offset a surface by a specified value."
    );
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("offset value");
    argList::validArgs.append("output surfaceFile");    
    argList::addBoolOption
    (
        "negative",
        "Offset in negative direction."
    );
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    scalar offset = args.argRead<scalar>(2);
    const fileName outFileName  = args[3];    
    const bool isNegative = args.optionFound("negative");
    if(isNegative)
    {
        offset *= -1;
    }

    Info<< "Reading surface from " << surfFileName << " ..." << endl << endl;

    triSurface surf(surfFileName);
    
    vectorField normals(surf.pointNormals());
    pointField points(surf.points());
    Map<label> map(surf.meshPointMap());
    
    forAll(points, pointI)
    {
        points[pointI] += normals[map[pointI]] * offset;
    }
    
    surf.movePoints(points);
    
    Info<< "Writing refined surface to " << outFileName << " ..." << endl;
    
    surf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

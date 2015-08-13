/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    modifyRegions

Description
    Modify regions created by splitMeshRegions -makeCellZones

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "Keyed.H"
#include "globalMeshData.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "pointZoneSet.H"



using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ZoneType>
void removeZone
(
    ZoneMesh<ZoneType, polyMesh>& zones,
    const word& setName
)
{
    label zoneID = zones.findZoneID(setName);

    if (zoneID != -1)
    {
        // Shuffle to last position
        labelList oldToNew(zones.size());
        label newI = 0;
        forAll(oldToNew, i)
        {
            if (i != zoneID)
            {
                oldToNew[i] = newI++;
            }
        }
        oldToNew[zoneID] = newI;
        zones.reorder(oldToNew);
        // Remove last element
        zones.setSize(zones.size()-1);
        zones.clearAddressing();
        zones.write();
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);
    #include "addDictOption.H"
    argList::addBoolOption
    (
        "merge",
        "Merge new regions into one single region."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedPolyMesh.H"

    const word dictName("modifyRegionsDict");
    #include "setSystemMeshDictionaryIO.H"
    IOdictionary modifyRegionsDict(dictIO);

    const bool relativeSize = readBool(modifyRegionsDict.lookup("relativeSize"));
    label minRegionSize;

    if (relativeSize)
    {
        const scalar ratio = readScalar(modifyRegionsDict.lookup("minRegionSize"));
        minRegionSize = int(mesh.nCells()*ratio + 0.5);
    }
    else
    {
        minRegionSize = readLabel(modifyRegionsDict.lookup("minRegionSize"));
    }

    const word defaultCellZone = modifyRegionsDict.lookup("defaultCellZone");
    const word cavitiesCellZone = modifyRegionsDict.lookup("cavitiesCellZone");
    const bool mergeRegions = readBool(modifyRegionsDict.lookup("mergeRegions"));
    const bool invertCavities = readBool(modifyRegionsDict.lookup("invertCavities"));
    const List<dictionary> cellZones = modifyRegionsDict.lookup("cellZones");

    Map<word> identifiers;
    forAll(cellZones, zoneI)
    {
        word zoneName = cellZones[zoneI].lookup("name");
        pointField insidePoints = cellZones[zoneI].lookup("insidePoints");
        forAll( insidePoints, pointI )
        {
            point insidePoint = insidePoints[pointI];
            label cellI = mesh.findCell(insidePoint);
            identifiers.insert(cellI, zoneName);
        }
    }
    const labelList insideCells = identifiers.toc();

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        const cellZoneMesh& cellZones = mesh.cellZones();
        const wordList cellZoneNames = cellZones.names();
        wordList solids;
        wordList regions;

        forAll( cellZoneNames, i )
        {
            word cellZoneName = cellZoneNames[i];
            if ( regExp("region.*").match(cellZoneName) )
                regions.append( cellZoneName );
            else
                solids.append( cellZoneName );
        }

        Info<< nl
            << "Remove multiply assigned cells from region cellSets:" 
            << nl ;

        forAll( regions, i )
        {
            word region = regions[i];
            autoPtr<topoSet> currentSet;
            const word setType = "cellSet";

            currentSet = topoSet::New
            (
                setType,
                mesh,
                region,
                IOobject::MUST_READ
            );

            Info<< "Read " << currentSet().type() << " "
                << region<< " with size "
                << returnReduce(currentSet().size(), sumOp<label>())
                << endl;

            forAll( solids, i )
            {
                word solid = solids[i];
                word sourceType = "cellToCell";
                dictionary sourceInfo;
                sourceInfo.add("set", solid);
                topoSetSource::setAction action = topoSetSource::DELETE;

                autoPtr<topoSetSource> source = topoSetSource::New
                (
                    sourceType,
                    mesh,
                    sourceInfo
                );

                source().applyToSet(action, currentSet());
                currentSet().write();
            }
        }

        const word setType = "cellSet";
        autoPtr<topoSet> cavitiesSet;
        cavitiesSet = topoSet::New
        (
            setType,
            mesh,
            cavitiesCellZone,
            10000
        );
        cavitiesSet().write();

        Info<< nl << "Assign regions to specified cellSets:" << nl;

        wordList setsToRemove;
        wordHashSet newSets;

        forAll( regions, i )
        {
            word region = regions[i];
            autoPtr<topoSet> regionSet;
            const word setType = "cellSet";
            regionSet = topoSet::New
            (
                setType,
                mesh,
                region,
                IOobject::MUST_READ
            );

            setsToRemove.append(region);

            label regionSize = regionSet().size();
            if (regionSize == 0)
            {
                continue;
            }

            word newCellSetName = cavitiesCellZone;
            if (regionSize > minRegionSize)
            {
                newCellSetName = defaultCellZone;
                newSets.insert(newCellSetName);
            }

            forAll( insideCells, i )
            {
                label insideCell = insideCells[i];
                if ( regionSet()[insideCell] )
                {
                    newCellSetName = identifiers[insideCell];
                    newSets.insert(newCellSetName);
                }
            }


            autoPtr<topoSet> newCellSet;
            newCellSet = topoSet::New
            (
                setType,
                mesh,
                newCellSetName,
                IOobject::READ_IF_PRESENT
            );

            Info<< "Create new cellSet " << newCellSetName << nl;

            word sourceType = "cellToCell";
            dictionary sourceInfo;
            sourceInfo.add("set", region);
            topoSetSource::setAction action = topoSetSource::ADD;
            autoPtr<topoSetSource> source = topoSetSource::New
            (
                sourceType,
                mesh,
                sourceInfo
            );

            source().applyToSet(action, newCellSet());
            newCellSet().write();

        }

        cavitiesSet = topoSet::New
        (
            setType,
            mesh,
            cavitiesCellZone,
            IOobject::MUST_READ
        );

        if (cavitiesSet().size() == 0)
        {
            setsToRemove.append(cavitiesCellZone);
        }
        else if (invertCavities)
        {
            Info<< "Invert set " << cavitiesCellZone << nl;
            cavitiesSet().invert(mesh.nCells());
            cavitiesSet().write();
        }

        Info<< nl;


        if (setsToRemove.size() != 0 )
        {
            IOobjectList objects
            (
                mesh,
                mesh.time().findInstance
                (
                    polyMesh::meshSubDir/"sets",
                    word::null,
                    IOobject::READ_IF_PRESENT,
                    mesh.facesInstance()
                ),
                polyMesh::meshSubDir/"sets"
            );

            forAll( setsToRemove, i )
            {
                word setName = setsToRemove[i];
                Info<< "Removing set " << setName << endl;
                if (objects.found(setName))
                {
                    fileName object = objects[setName]->objectPath();
                    rm(object);
                }
                removeZone
                (
                    const_cast<cellZoneMesh&>(mesh.cellZones()),
                    setName
                );
            }
        }


        wordList newSetsList = newSets.toc();
        forAll( newSetsList, i )
        {
            const word newSet = newSetsList[i];
            const word setType = "cellZoneSet";
            autoPtr<topoSet> newCellZoneSet;
            newCellZoneSet = topoSet::New
            (
                setType,
                mesh,
                newSet,
                10000
            );

            Info<< "Create new cellZoneSet " << newSet << nl;

            word sourceType = "setToCellZone";
            dictionary sourceInfo;
            sourceInfo.add("set", newSet);
            topoSetSource::setAction action = topoSetSource::ADD;
            autoPtr<topoSetSource> source = topoSetSource::New
            (
                sourceType,
                mesh,
                sourceInfo
            );

            source().applyToSet(action, newCellZoneSet());
            newCellZoneSet().write();
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

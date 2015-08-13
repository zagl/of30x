
fvMesh& mesh = solidRegions[i];
const dictionary& subcycleControl = mesh.solutionDict().subDict("subcycleControl");

double maxDeltaTForSolidRegion =
    subcycleControl.lookupOrDefault<double>("maxDeltaTForSolidRegion", 1.0);

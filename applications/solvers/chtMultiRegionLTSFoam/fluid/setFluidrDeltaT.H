{
    const dictionary& pimpleDict = runTime.controlDict();

    scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );
    
    forAll(fluidRegions, i)
    {
        fvMesh& mesh = fluidRegions[i];
    
        volScalarField rDeltaT0("rDeltaT0", rDeltaTFluid[i]);

        // Set the reciprocal time-step from the local Courant number
        rDeltaTFluid[i].dimensionedInternalField() = max
        (
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            fvc::surfaceSum(mag(phiFluid[i]))().dimensionedInternalField()
           /((2*maxCo)*mesh.V()*rhoFluid[i].dimensionedInternalField())
        );

        // Update tho boundary values of the reciprocal time-step
        rDeltaTFluid[i].correctBoundaryConditions();

        Info<< "Flow time scale min/max = "
            << gMin(1/rDeltaTFluid[i].internalField())
            << ", " << gMax(1/rDeltaTFluid[i].internalField()) << endl;

        if (rDeltaTSmoothingCoeff < 1.0)
        {
            fvc::smooth(rDeltaTFluid[i], rDeltaTSmoothingCoeff);
        }

        Info<< "Smoothed flow time scale min/max = "
            << gMin(1/rDeltaTFluid[i].internalField())
            << ", " << gMax(1/rDeltaTFluid[i].internalField()) << endl;

        // Limit rate of change of time scale
        // - reduce as much as required
        // - only increase at a fraction of old time scale
        if
        (
            rDeltaTDampingCoeff < 1.0
         && runTime.timeIndex() > runTime.startTimeIndex() + 1
        )
        {
            rDeltaTFluid[i] =
                rDeltaT0
               *max(rDeltaTFluid[i]/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

            Info<< "Damped flow time scale min/max = "
                << gMin(1/rDeltaTFluid[i].internalField())
                << ", " << gMax(1/rDeltaTFluid[i].internalField()) << endl;
        }
    }
}

{
    const dictionary& pimpleDict = runTime.controlDict();

    scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

    scalar maxDi
    (
        pimpleDict.lookupOrDefault<scalar>("maxDi", 10.0)
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

    forAll(solidRegions, i)
    {
        fvMesh& mesh = solidRegions[i];
        const solidThermo& thermo = thermos[i];

        tmp<volScalarField> magKappa;
        if (thermo.isotropic())
        {
            magKappa = thermo.kappa();
        }
        else
        {
            magKappa = mag(thermo.Kappa());
        }
        const volScalarField& kappa = magKappa();

        tmp<volScalarField> tcp = thermo.Cp();
        const volScalarField& cp = tcp();

        tmp<volScalarField> trho = thermo.rho();
        const volScalarField& rho = trho();

        volScalarField rDeltaT0("rDeltaT0", rDeltaTSolid[i]);

        surfaceScalarField kapparhoCpbyDelta
        (
            mesh.surfaceInterpolation::deltaCoeffs()
          * mesh.surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(kappa)
          / fvc::interpolate(cp*rho)
          / maxDi
        );

        // Set the reciprocal time-step from the local Courant number
        rDeltaTSolid[i].dimensionedInternalField() = max
        (
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            fvc::surfaceSum(kapparhoCpbyDelta)().dimensionedInternalField()
        );

        // Update tho boundary values of the reciprocal time-step
        rDeltaTSolid[i].correctBoundaryConditions();

        Info<< "Flow time scale min/max = "
            << gMin(1/rDeltaTSolid[i].internalField())
            << ", " << gMax(1/rDeltaTSolid[i].internalField()) << endl;

        if (rDeltaTSmoothingCoeff < 1.0)
        {
            fvc::smooth(rDeltaTSolid[i], rDeltaTSmoothingCoeff);
        }

        Info<< "Smoothed flow time scale min/max = "
            << gMin(1/rDeltaTSolid[i].internalField())
            << ", " << gMax(1/rDeltaTSolid[i].internalField()) << endl;

        // Limit rate of change of time scale
        // - reduce as much as required
        // - only increase at a fraction of old time scale
        if
        (
            rDeltaTDampingCoeff < 1.0
         && runTime.timeIndex() > runTime.startTimeIndex() + 1
        )
        {
            rDeltaTSolid[i] =
                rDeltaT0
               *max(rDeltaTSolid[i]/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

            Info<< "Damped flow time scale min/max = "
                << gMin(1/rDeltaTSolid[i].internalField())
                << ", " << gMax(1/rDeltaTSolid[i].internalField()) << endl;
        }
    }
}

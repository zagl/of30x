
    PtrList<volScalarField> rDeltaTFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        
        rDeltaTFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rDeltaT",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        
    }

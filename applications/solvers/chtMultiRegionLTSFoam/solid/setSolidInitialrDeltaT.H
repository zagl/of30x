
    PtrList<volScalarField> rDeltaTSolid(solidRegions.size());

    forAll(solidRegions, i)
    {
        
        rDeltaTSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rDeltaT",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        
    }

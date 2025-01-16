// Initialise fluid field pointer lists
PtrList<fluidReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> pRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::momentumTransportModel>
    turbulenceFluid(fluidRegions.size());
PtrList<fluidReactionThermophysicalTransportModel>
    thermophysicalTransportFluid(fluidRegions.size());
PtrList<combustionModel> reactionFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<Foam::fvModels> fvModelsFluid(fluidRegions.size());
PtrList<fvConstraints> fvConstraintsFluid(fluidRegions.size());

PtrList<dimensionedScalar> epsReactor(fluidRegions.size());
PtrList<dimensionedScalar> KcReactor(fluidRegions.size());
PtrList<dimensionedScalar> rhoCpcReactor(fluidRegions.size());
PtrList<dimensionedScalar> dpReactor(fluidRegions.size());
PtrList<dimensionedScalar> AoVReactor(fluidRegions.size());

PtrList<dimensionedScalar> a1Reactor(fluidRegions.size());
PtrList<dimensionedScalar> a2Reactor(fluidRegions.size());
PtrList<dimensionedScalar> a3Reactor(fluidRegions.size());

PtrList<volScalarField> ReReactor(fluidRegions.size());
PtrList<volScalarField> PrReactor(fluidRegions.size());
PtrList<volScalarField> NuReactor(fluidRegions.size());
PtrList<volScalarField> hcReactor(fluidRegions.size());

PtrList<volScalarField> TcReactor(fluidRegions.size());
PtrList<volScalarField> zFReactor(fluidRegions.size());

List<bool> activeFluid(fluidRegions.size(), false);
List<bool> activeCatFluid(fluidRegions.size(), false);
PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, fluidReactionThermo::New(fluidRegions[i]).ptr());

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimLength, 0)
        )
    );

    Info<< "    Adding to pRefFluid\n" << endl;
    pRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "pRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimPressure, 0)
        )
    );

    dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::momentumTransportModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to thermophysicalTransport\n" << endl;
    thermophysicalTransportFluid.set
    (
        i,
        fluidReactionThermophysicalTransportModel::New
        (
            turbulenceFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        combustionModel::New(thermoFluid[i], turbulenceFluid[i])
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();


    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            fluidRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    transportProperties.readIfPresent("inertMixture", activeFluid[i]);
    transportProperties.readIfPresent("CatalyticActive", activeCatFluid[i]);

    Info<< "Reading porousity eps\n" << endl;
    epsReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("eps")
	    )
    );

    Info<< "Reading catalyst conductivity\n" << endl;
    KcReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("Kc")
	    )
    );

    Info<< "Reading catalyst heat capacity\n" << endl;
    rhoCpcReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("rhoCpc")
	    )
    );

    Info<< "Reading pellet diameter\n" << endl;
    dpReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("dp")
	    )
    );

    Info<< "Reading pellet diameter\n" << endl;
    a1Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a1")
	    )
    );

    Info<< "Reading pellet diameter\n" << endl;
    a2Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a2")
	    )
    );

    Info<< "Reading pellet diameter\n" << endl;
    a3Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a3")
	    )
    );

    Info<< "Reading pellet diameter\n" << endl;
    AoVReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("AoV")
	    )
    );


    Info<< "    Adding to ReReactor\n" << endl;
    ReReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Re",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i],
	    dimensionedScalar("Re", dimless, scalar(0.0))
        )
    );

    Info<< "    Adding to PrReactor\n" << endl;
    PrReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Pr",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i],
	    dimensionedScalar("Pr", dimless, scalar(0.7))
        )
    );

    Info<< "    Adding to NuReactor\n" << endl;
    NuReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Nu",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i],
	    dimensionedScalar("Nu", dimless, scalar(0.0))
        )
    );

    Info<< "    Adding to ReReactor\n" << endl;
    hcReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "hc",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i],
	    dimensionedScalar("hc", dimensionSet(1, 0, -3, -1, 0, 0, 0), scalar(0.0))
        )
    );

    Info<< "    Adding to TcReactor\n" << endl;
    TcReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Tc",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i]
        )
    );

    Info<< "    Adding to zFReactor\n" << endl;
    zFReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "zF",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    fluidRegions[i],
	    dimensionedScalar("zF", dimless, scalar(0.0))
        )
    );

    rhoMax.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMax",
                transportProperties,
                dimDensity,
                GREAT
            )
        )
    );

    rhoMin.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMin",
                transportProperties,
                dimDensity,
                0
            )
        )
    );    

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvModelsFluid\n" << endl;
    fvModelsFluid.set
    (
        i,
        new Foam::fvModels(fluidRegions[i])
    );

    Info<< "    Adding fvConstraintsFluid\n" << endl;
    fvConstraintsFluid.set
    (
        i,
        new fvConstraints(fluidRegions[i])
    );

    turbulenceFluid[i].validate();
}

    const fvMesh& mesh = fluidRegions[i];

    combustionModel& reaction = reactionFluid[i];

    fluidReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::momentumTransportModel& turbulence = turbulenceFluid[i];
    fluidReactionThermophysicalTransportModel& thermophysicalTransport =
        thermophysicalTransportFluid[i];

    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const dimensionedScalar& pRef = pRefFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    Foam::fvModels& fvModels = fvModelsFluid[i];
    Foam::fvConstraints& fvConstraints = fvConstraintsFluid[i];

    #include "checkRadiationModel.H"

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimple = pimples.pimple(i);

    pressureReference& pressureReference = pressureReferenceFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];

    // This solver does not support moving mesh but it uses the pressure
    // equation of one which does, so we need a dummy face-momentum field
    autoPtr<surfaceVectorField> rhoUf(nullptr);

    volScalarField& Tc=TcReactor[i];
    volScalarField& zF=zFReactor[i];
    volScalarField& Re=ReReactor[i];
    volScalarField& Nu=NuReactor[i];
    volScalarField& hc=hcReactor[i];
    volScalarField& Pr=PrReactor[i];

    const dimensionedScalar& dp = dpReactor[i]; 
    const dimensionedScalar& eps = epsReactor[i];
    const dimensionedScalar& Kc = KcReactor[i];
    const dimensionedScalar& rhoCpc = rhoCpcReactor[i];
    const dimensionedScalar& AoV = AoVReactor[i];

    const dimensionedScalar& a1 = a1Reactor[i];
    const dimensionedScalar& a2 = a2Reactor[i];
    const dimensionedScalar& a3 = a3Reactor[i]; 
    
    const bool active = activeFluid[i];
    const bool activeCat = activeCatFluid[i];
    if (activeCat)
    {
        word cellSetName = "catalyst";

        label zoneID = mesh.cellZones().findZoneID(cellSetName);

        if (zoneID == -1)
        {
        FatalErrorIn("yourFunctionName")
            << "Cannot find cellZone " << cellSetName << endl
            << "Valid cellZones are " << mesh.cellZones().names()
            << exit(FatalError);
        }

        const labelList& cells = mesh.cellZones()[zoneID];

        Info << "Cells in cellzone " << cellSetName << ":" << endl;
        forAll(cells, i)
        {
        const label cell = cells[i];
        zF[cell]=1;
        }
    }

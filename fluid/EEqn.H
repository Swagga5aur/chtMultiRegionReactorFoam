{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + (
            mvConvection.valid()
          ? mvConvection->fvmDiv(phi, he)
          : fvm::div(phi, he)
        )
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection.valid()
            ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
            : fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport.divq(he)
     ==
        rho*(U&g)
      + reaction.Qdot()*(1-zF)
      + fvModels.source(rho, he)
      - fvc::Sp(hc*AoV,thermo.T())*zF
      + hc*AoV*Tc*zF
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
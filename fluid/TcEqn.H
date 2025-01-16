if (activeCat)
{
    Re = dp*mag(U)*rho/thermo.mu();
    Pr = thermo.mu()*thermo.Cp()/thermo.kappa();
    Nu = a1/eps*pow(Pr,a3)*pow(Re,a2); //a1*pow(Re,a2)*pow(Pr,a3);	a1/eps*pow(Pr,a3)*pow(Re,a2);
    hc = Nu*thermo.kappa()/dp;	//Nu*thermo.kappa()/dp

    fvScalarMatrix TcEqn
    (
	(1-eps)*rhoCpc*fvm::ddt(Tc)
	-fvm::laplacian((scalar(1)-eps)*Kc,Tc,"laplacian(DT,T)")
	==
	 reaction.Qdot()
	-fvm::Sp(hc*AoV,Tc) 
	+hc*AoV*thermo.T()

    );
	
    TcEqn.relax();

    TcEqn.solve();

    Tc = Tc*zF + thermo.T()*(1-zF);

    Info<< "Min/max Tcatalyst:" << min(Tc).value() << ' '
        << max(Tc).value() << endl;
} else
{
    Tc = thermo.T();
}
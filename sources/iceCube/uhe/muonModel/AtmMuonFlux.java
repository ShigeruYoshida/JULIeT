package iceCube.uhe.muonModel;

import java.io.*;
import java.util.*;

import numRecipes.*;

/**
<pre>
UHE Atmospheric Muon fluxes based on the following model
are calculated from the table.

</pre>
<UL>
  <DT> <cite>Atmospheric Neutrino Flux above 1GeV</cite> 
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> V. Agrawal et al ,</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://publish.aps.org/abstract/PRD/v53/p1314">
            Phys.Rev.Rev. <b>D 55</b> 1314 (1996)</a></TD>
       </TR>
       </TABLE><BR>
</UL>

*/


public class AtmMuonFlux extends ParticleFlux {

    private int numberOfFlavor = 2;
    private static final double ln10 = Math.log(10.0);
    private int dataNumber = 1; // Initial value. Should be changed later.

    private double[] EFluxArray = new double [6];
    private double[] EFluxArray5TeV = new double [6];
    private double[] cosThetaArray = new double [6];


    /** Constructor: Energy Flux dF/dE E^2 [GeV/cm^2 sec sr] at 5 TeV is given.
    */
    public AtmMuonFlux() {

	cosThetaArray[0] = 1.0;
	EFluxArray5TeV[0] = 9.80e-6;
	cosThetaArray[1] = 0.75;
	EFluxArray5TeV[1] = 1.445e-5;
	cosThetaArray[2] = 0.50;
	EFluxArray5TeV[2] = 1.92e-5;
	cosThetaArray[3] = 0.25;
	EFluxArray5TeV[3] = 3.66e-5;
	cosThetaArray[4] = 0.15;
	EFluxArray5TeV[4] = 3.939e-5;
	cosThetaArray[5] = 0.05;
	EFluxArray5TeV[5] = 5.5e-5;

    }


    /** 
	calculate the differential Energy Flux [GeV /cm^2 sec sr] 

	<pre>
	logEnergy [GeV]

	cos(Zenith Angle)
	</pre>

	E^2dF/dE ~ E^(-1.7) is assumed for extrapolation to UHE.

    */
    public double getEFlux(double logEnergy, double cosTheta){

	double energyFactor = Math.pow(10.0,logEnergy)/5.0e3; // Normalized to 5TeV=5000GeV;
	for(int icos=0;icos<6;icos++){
	    EFluxArray[icos] = EFluxArray5TeV[icos]*Math.pow(energyFactor,-1.7);
	}


	double EFlux = 
	    Interpolation.mThPolynominalInterpolate(cosThetaArray,
	    EFluxArray,6,cosTheta,2);

	
	return(EFlux);
    }


    /** 
	<pre>
	calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 

	logEnergy [GeV]

	cos(Zenith Angle)
	</pre>

    */
    public double getDFDLogE(double logEnergy, double cos_theta){

	double EFlux = getEFlux(logEnergy,cos_theta);

	double energy = Math.pow(10.0,logEnergy); // [GeV]
	double flux = EFlux/energy*ln10;

	return(flux);
    }

}


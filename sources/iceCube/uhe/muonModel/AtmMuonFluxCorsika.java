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


public class AtmMuonFluxCorsika extends ParticleFlux {

    private int numberOfFlavor = 2;
    private static final double ln10 = Math.log(10.0);
    private int dataNumber = 1; // Initial value. Should be changed later.

    private double[] EFluxArray = new double [6];
    private double[] EFluxArray10PeV = new double [6];
    private double[] cosThetaArray = new double [6];


    /** Constructor: Energy Flux dF/dE E^2 [GeV/cm^2 sec sr] at 10^7 GeV is given. <---

	Also, probably replace EFluxArray10PeV with EFluxArray10PeV
    */
    public AtmMuonFluxCorsika() {

	cosThetaArray[0] = 1.0;
	cosThetaArray[1] = 0.75;
	cosThetaArray[2] = 0.50;
	cosThetaArray[3] = 0.25;
	cosThetaArray[4] = 0.15;
	cosThetaArray[5] = 0.05;

	double fixedEnergy=1.e7;
	for(int i=0;i<=5;i++){
	    double costh=cosThetaArray[i];

	    /*
	      costh correction needs to go in here, e.g.
	    */

	    double h=73, R=6.4e3, xx;
	    xx=(R-h)*costh;
	    costh=h/(-xx+Math.sqrt(xx*xx-h*h+2*R*h));

	    EFluxArray10PeV[i] = 0.14*55.096*Math.pow(fixedEnergy,(-3.002+2.0))*
		(1.0/(1.0+1.1*fixedEnergy*(costh)/115.0)+0.054/(1.0+1.1*fixedEnergy*(costh)/850.0));
	}

    }


    /** 
	calculate the differential Energy Flux [GeV /cm^2 sec sr] 

	<pre>
	logEnergy [GeV]

	cos(Zenith Angle)
	</pre>

	E^2dF/dE ~ E^(-2.002) is assumed for extrapolation to UHE. <---

    */
    public double getEFlux(double logEnergy, double cosTheta){

	double energyFactor = Math.pow(10.0,logEnergy)/1.0e7; // Normalized to 10^7GeV;

	for(int icos=0;icos<6;icos++){
	    EFluxArray[icos] = EFluxArray10PeV[icos]*Math.pow(energyFactor,-2.002);
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



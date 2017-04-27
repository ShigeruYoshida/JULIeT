package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** 
Calculate the neutrino detection effective area [km^2 sr]  
(neutrino interaction probability convoluted)
by running PropagationMatrixFlux with EffAreaTable class.
Written by S.Yoshida 2007 April 19th
*/

public class DrawNuEffAreaByAreaTable {


    public static void main(String[] args) throws IOException{


	double solidAngle = 4.0*Math.PI; 
	PropagationMatrixFlux iceFlux = new PropagationMatrixFlux();


	// Interactive session to set the inIceParticle flavor and doublets
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 

	boolean interactionEnd = false;
	do{
	    System.err.print("inIceParticle flavor ->"); 
	    buffer   = d.readLine(); 
	    int flavor = Integer.valueOf(buffer).intValue();

	    System.err.print("inIceParticle doublet ->"); 
	    buffer   = d.readLine(); 
	    int doublet = Integer.valueOf(buffer).intValue();

	    iceFlux.addInIceParticle(flavor,doublet);

	    System.err.print("No more inIce particles [yes(1)/no(0)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) interactionEnd = true; 
	}while(!interactionEnd);

	iceFlux.calculateYield();

	    
        System.out.println("titx Energy [GeV]"); 
        System.out.println("tity Area [km^2!])"); 
        System.out.println("scal 1.0e6 1.0e11 3.0e-7 0.1"); 
 

	double epsilon = 1.0e-4; // round-off margin for binning
	double logNeutrinoEnergy = Particle.getLogEnergyMinimum() + epsilon;
	double logNeutrinoMaxEnergy = Particle.getLogEnergyMinimum()
	+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	while(logNeutrinoEnergy < logNeutrinoMaxEnergy){ // E < 10^12 GeV

            double area = iceFlux.getYield(logNeutrinoEnergy)*3.0/ 
                (iceFlux.observationTime*1.0e10*solidAngle); 
            double nuEnergy = Math.pow(10.0,(logNeutrinoEnergy-epsilon)); 
            System.out.println("data " + nuEnergy + " 0.0 " + 
                               area  + " 0.0"); 
	    logNeutrinoEnergy += Particle.getDeltaLogEnergy();


        } 
 
        System.out.println("logx"); 
        System.out.println("logy"); 
        System.out.println("join"); 
        System.out.println("disp"); 
        System.out.println("endg"); 

    }


}

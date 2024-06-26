package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of Leptonic Glashow Resonance Interactions */
public class MakeGlashowResonanceLeptonicMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor = 0;         // Incident particle flavor,  e->0
        int doublet =0;         // Incident particle doublet, neutrinos->0
        int producedFlavor = 0; // Produced particle flavor,  e->0, mu->1, tau->2
        int material = 0;       // Ice->0, Rock->1
        Boolean is_nlo = false; // Whether to use the LO or NLO Glashow xsection
        double energy = 1.0e9;  // 1EeV= 10^9 GeV
        double epsilon = 1.0e-6;
        double energyCut = 1.0e4; // 10000 GeV


        if(args.length!=4){
            System.out.println("Usage: MakeGlashowResonanceLeptonicMtx file-name producedFlavor material is-nlo");
        System.exit(0);
        }else{
            fileName = args[0];
            producedFlavor = Integer.valueOf(args[1]).intValue();
            material = Integer.valueOf(args[2]).intValue();
            is_nlo = Boolean.valueOf(args[3]);
        }

    numRecipes.Integration.setRelativeAccuracy(epsilon);
    // integration accuracy for save CPUtime

    // Generate the ParticlePoint class.
    ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

    for(int i=0;i<s.NumberOfSpecies[material];i++){
        System.out.println("Charge " + s.getCharge(i));
        System.out.println("Atomic Number " + s.getAtomicNumber(i));
    }

    //Generate object of Leptonic Glashow Resonance interaction.
    if(is_nlo){
        GlashowResonanceLeptonicNLO grLeptonNLO = new GlashowResonanceLeptonicNLO(s, producedFlavor);
        System.err.println(grLeptonNLO.interactionName());

        //Generate object of the Interaction Matrix
        InteractionsMatrix grLeptonMtx = new InteractionsMatrix(grLeptonNLO);
        int iLogE;
        for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
            grLeptonNLO.setIncidentParticleEnergy(iLogE);
            if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

            // Total Cross Section
            grLeptonMtx.setSigmaMatrix(iLogE);
            System.err.println("  Total cross section done");

            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<Particle.getDimensionOfLogEnergyMatrix();jLogE++){
                grLeptonMtx.setTransferMatrix(iLogE,jLogE);
            }

        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(grLeptonMtx, out);
        out.close();
    }
    else{
        GlashowResonanceLeptonic grLepton = new GlashowResonanceLeptonic(s, producedFlavor);
        System.err.println(grLepton.interactionName());

        //Generate object of the Interaction Matrix
        InteractionsMatrix grLeptonMtx = new InteractionsMatrix(grLepton);
        int iLogE;
        for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
            grLepton.setIncidentParticleEnergy(iLogE);
            if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

            // Total Cross Section
            grLeptonMtx.setSigmaMatrix(iLogE);
            System.err.println("  Total cross section done");


            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<Particle.getDimensionOfLogEnergyMatrix();jLogE++){
                grLeptonMtx.setTransferMatrix(iLogE,jLogE);
            }
        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(grLeptonMtx, out);
        out.close();
    }
}
}

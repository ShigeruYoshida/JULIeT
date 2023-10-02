package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of Hadronic Glashow Resonance Interactions. */
public class MakeGlashowResonanceHadronicMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor = 0;         // Incident particle flavor,  e->0
        int doublet =0;         // Incident particle doublet, neutrinos->0
        int material = 0;       // Ice->0, Rock->1
        Boolean is_nlo = false; // Whether to use the LO or NLO Glashow xsection
        double energy = 1.0e9;  // 1EeV= 10^9 GeV
        double epsilon = 1.0e-3;
        double energyCut = 1.0e4; // 10000 GeV


        if(args.length!=3){
            System.out.println("Usage: MakeGlashowResonanceHadronicMtx file-name material is-nlo");
        System.exit(0);
        }else{
            fileName = args[0];
            material = Integer.valueOf(args[1]).intValue();
            is_nlo = Boolean.valueOf(args[2]);
        }


    // Generate the ParticlePoint class.
    ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

    for(int i=0;i<s.NumberOfSpecies[material];i++){
        System.out.println("Charge " + s.getCharge(i));
        System.out.println("Atomic Number " + s.getAtomicNumber(i));
    }

    numRecipes.Integration.setRelativeAccuracy(epsilon);
    // integration accuracy for save CPUtime

    //Generate object of Hadronic Glashow Resonance interaction.
    if(is_nlo){
        GlashowResonanceHadronicNLO grHadronNLO = new GlashowResonanceHadronicNLO(s);
        System.err.println(grHadronNLO.interactionName( ));

        //Generate object of the Interaction Matrix
        GlashowResonanceHadronicMatrix grHadronMtx = new GlashowResonanceHadronicMatrix(grHadronNLO);

        int iLogE;
        for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
            grHadronNLO.setIncidentParticleEnergy(iLogE);
            if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

            // Total Cross Section
            grHadronMtx.setSigmaMatrix(iLogE);
            System.err.println("  Total cross section done");

            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<Particle.getDimensionOfLogEnergyMatrix();jLogE++){
                grHadronMtx.setTransferMatrix(iLogE,jLogE);
            }

        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(grHadronMtx, out);
        out.close();
    }
    else{
        GlashowResonanceHadronic grHadron = new GlashowResonanceHadronic(s);
        System.err.println(grHadron.interactionName( ));

        //Generate object of the Interaction Matrix
        GlashowResonanceHadronicMatrix grHadronMtx = new GlashowResonanceHadronicMatrix(grHadron);

        int iLogE;
        for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
            grHadron.setIncidentParticleEnergy(iLogE);
            if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

            // Total Cross Section
            grHadronMtx.setSigmaMatrix(iLogE);
            System.err.println("  Total cross section done");

            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<Particle.getDimensionOfLogEnergyMatrix();jLogE++){
            grHadronMtx.setTransferMatrix(iLogE,jLogE);
            //if(jLogE%100==0) System.err.println("  Transfer Matrix " + jLogE);
            }

        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(grHadronMtx, out);
        out.close();
    }
}
}

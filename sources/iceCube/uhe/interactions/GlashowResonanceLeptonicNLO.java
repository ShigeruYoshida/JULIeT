package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;
import numRecipes.Interpolation;

/** 
  The Glashow Resonance reaction with W into the leptonic decay
   <pre>
   \bar{nu_e} + e^{-1} -> \bar{nu_l} + l ^{-1}
   </pre>
  The approximation that
   <pre>
      1. Masses of the produced leptons are negligible
      2. Energy of incoming anti electron neutrino is by far higher
         than the electron mass.
   </pre>
   has been introduced to calculate the differential cross section

   The inelasiticity parameter y is here difined as
   <pre>
     y = 1 - E_{l^{-1}}/E_{\bar{\nu_e}}
   </pre>

*/

class DataFileLoader {
    public static CrossSectionData loadDataFile(String filePath) {
        List<Double> energy = new ArrayList<>();
        List<Double> crossSection = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;

            while ((line = reader.readLine()) != null) {
                String[] columns = line.split("\\s+");

                if (columns.length == 2) {
                    double energyValue = Double.parseDouble(columns[0]);
                    double crossSectionValue = Double.parseDouble(columns[1]);

                    energy.add(energyValue);
                    crossSection.add(crossSectionValue);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[] energyArray = listToDoubleArray(energy);
        double[] crossSectionArray = listToDoubleArray(crossSection);

        return new CrossSectionData(energyArray, crossSectionArray);
    }

    private static double[] listToDoubleArray(List<Double> list) {
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
}

class CrossSectionData {
        private double[] energy;
        private double[] crossSection;

        public CrossSectionData(double[] energy, double[] crossSection) {
            this.energy = energy;
            this.crossSection = crossSection;
        }

        public double[] getEnergy() {
            return energy;
        }

        public double[] getCrossSection() {
            return crossSection;
        }
}


public class  GlashowResonanceLeptonicNLO extends Interactions implements Function{
    private String filePath = "sources/iceCube/uhe/interactions/xsec_GR_full.dat";

    double[] energy_data;
    double[] crossSection_data;

    protected boolean isPerNucleon = true; // Cross secion is given as per target nucleon
    protected double chargePerNucleon = 0.0;  // number of electrons per nucleon
                                              // calculated in the constructor.
    protected final double fractionOfLeptonicDecay = 0.1086; // from PDG 2023

    /** Constructor: Register the ParticlePoint classes 
         and the prodocued flavor - 0 for e, 1 for mu, 2 for tau
    */
    public GlashowResonanceLeptonicNLO(ParticlePoint s, int flavor){
        super(new Particle(0,0), s, flavor); 
        // Progagating particle is an electron neutrino

        CrossSectionData data = DataFileLoader.loadDataFile(filePath);
        energy_data = data.getEnergy();
        crossSection_data = data.getCrossSection();


        // calculate the number of electrons per nucleon
        double massNumber = 0.0;
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
            chargePerNucleon += (double )(s.getNumberOfAtoms(i))*s.getCharge(i);
            massNumber += (double )(s.getNumberOfAtoms(i))*s.getAtomicNumber(i);
        }
        chargePerNucleon = chargePerNucleon/massNumber;

    }

    /** Checking the particle kind involved with
    a given interaction. Only electron neutrino
    is allowed to be involved with the Glashow resonance
    in the medium. This is the overriden method from Interactions.java
    */
    public boolean isValidParticle(Particle p){
        if(p.getDoublet()==0 && p.getFlavor()==0) return true;
        else return false;
    }

    /**
       Calculate the differential cross section as the one per nucleon.
       for instance in the case of ice, the number of electrons per nuclei
       is approximately (1 x 8 + 2 x1) - corresponding
       to H_2 0 !! -- is multiplied to the cross section.
    */
    public void calculateCrossSectionAsPerNucleon(){
        isPerNucleon = true;
    }

    /**
       Calculate the differential cross section as the one per electron.
       This is simple - because the Glashow resonance is subject to
       electron.
    */
    public void calculateCrossSectionAsPerElectron(){
        isPerNucleon = false;
    }

    /**
       return dSigma/dy [cm^2] where y = 1 - - E_{l^{-1}}/E_{\bar{\nu_e}}
    */
    public double getDSigmaDy(double y){
        if(!isValidInelasticity(y)) return 0.0;
        return 3*getSigma()*y*y;
    }

    /** return total cross section [cm^2] */
    public double getSigma(){
        if (energy < energy_data[0] || energy > energy_data[energy_data.length-1]){
            return 0;
        }
        double sigma = Interpolation.mThPolynominalInterpolate(
            energy_data, crossSection_data, energy_data.length, energy, 4);
        if(!isPerNucleon){
            return sigma * fractionOfLeptonicDecay;
        }
        else
            return chargePerNucleon * sigma * fractionOfLeptonicDecay;
    }

    /** Checking the range of the given inelasticity y 
    that is determined in an individual interaction channel.
    */
    public boolean isValidInelasticity(double y){ 
        if((getYmin()+roundOffError)<= y && 
           y <= getYmax()) return true;
        else return false;
    }

    public double getYmin(){
        return 0.0;
    }

    public double getYmax(){
        return 1.0;
    }

    public String interactionName(){
        String channel = "Glashow Resonance Leptonic ";
        return channel;
    }
}

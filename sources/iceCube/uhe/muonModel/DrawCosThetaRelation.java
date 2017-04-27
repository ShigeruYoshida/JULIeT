package iceCube.uhe.muonModel;

import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawCosThetaRelation {


   public static void main(String[] args) throws IOException {

       AtmMuonBundleFlux muonFlux = new AtmMuonBundleFlux();

       System.out.println("titx Cos(Zenith at Ice3)");
       System.out.println("tity Cos(Zenith at the surface)");
       System.out.println("gwin 0.2 0.9 0.2 0.9");
       System.out.println("scal -1.0 1.0 -1.0 1.0");

       double cosTheta = -1.0;
       while(cosTheta<=1.0){
	   double cosThetaSurface = 
	       muonFlux.getCosineOfZenithAngleAtEarthSurface(cosTheta);
	   System.out.println("data " + cosTheta + " 0.0 " + cosThetaSurface + " 0.0");
	   cosTheta += 0.05;
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");
	System.out.println("endg");
	    
   }
}

package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

public class TableAdjust {

    static double[ ][ ] yDsigmaArray = new double[32][14];
    static double[ ] logEArray = new double[32];
    static double[ ] logDyArray = new double[14];

    public static void main(String[] args) throws IOException {

	int iLogE;
	String inFile = null;
	String outFile = null;

	if(args.length!=2){
            System.out.println("Usage: TableAdjust file-name(in) file-name(out)");
	    System.exit(0);
        }else{
            inFile = args[0];
            outFile = args[1];
        }

         
        DataInputStream in =
            new DataInputStream(new FileInputStream(inFile));


        int jLogY;
        for(jLogY=0;jLogY<14;jLogY++) logDyArray[jLogY] =in.readDouble();

        for(iLogE=0;iLogE<=31;iLogE++){
            logEArray[iLogE] = in.readDouble( );
            for(jLogY=0;jLogY<14;jLogY++){
                double ySigma = in.readDouble();
                if(ySigma != Double.NEGATIVE_INFINITY && ySigma != Double.NaN){
                    yDsigmaArray[iLogE][jLogY] = ySigma;
                }else{
                    yDsigmaArray[iLogE][jLogY] = yDsigmaArray[iLogE][jLogY-1];
                }
            }
        }

	in.close( );

	for(iLogE=0;iLogE<=31;iLogE++){
	    for(jLogY=1;jLogY<14;jLogY++){
              if(yDsigmaArray[iLogE][jLogY] <= -40.0){
                  int kLogY;int offset = 0;
		  do{
		      kLogY = replace(iLogE,jLogY,iLogE-offset);
		      offset++;
		  }while(kLogY==0);
              }
	    }

	}


	DataOutputStream out =
	    new DataOutputStream(new FileOutputStream(outFile));

	for(jLogY=0;jLogY<14;jLogY++) out.writeDouble(logDyArray[jLogY]);
	for(iLogE=0;iLogE<32;iLogE+=5){
	    out.writeDouble(logEArray[iLogE]);
	    for(jLogY=0;jLogY<14;jLogY++){
		out.writeDouble(yDsigmaArray[iLogE][jLogY]);
	    }
	}

	out.writeDouble(logEArray[31]);
	for(jLogY=0;jLogY<14;jLogY++){
	    out.writeDouble(yDsigmaArray[31][jLogY]);
	}
	out.close( );
    }



    public static int replace(int iLogE,int jLogY, int kLogE){

	int kLogY;
	for(kLogY=jLogY;kLogY<14;kLogY++){
	    if(yDsigmaArray[kLogE][kLogY] > -40.0){
		System.err.println("iLogE " + iLogE + " jLogY " + jLogY +
				   " kLogE " + kLogE + " kLogY " + kLogY);
		yDsigmaArray[iLogE][jLogY] = yDsigmaArray[kLogE][kLogY];
		break;
	    }
	}
	if(kLogY == 14){
	    for(kLogY=jLogY;kLogY>0;kLogY--){
		if(yDsigmaArray[kLogE][kLogY] > -40.0){
		    System.err.println("iLogE " + iLogE + " jLogY " + jLogY +
				   " kLogE " + kLogE + " kLogY " + kLogY);
		    yDsigmaArray[iLogE][jLogY] = yDsigmaArray[kLogE][kLogY];
		    break;
		}
	    }
	}


	return kLogY;
    }

}



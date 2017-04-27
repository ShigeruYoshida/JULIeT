package numRecipes;

import java.io.*;

/** 
   <pre>
   Ported to Java by Shigeru Yoshida for the IceCube MC.
                     2002.12.1                          
                     syoshida@hepburn.s.chiba-u.ac.jp
   Before using, initialize the state by the constrctor RandomDouble(long)
   in this Java version.

   Blow is the notice of the original C source code.
   -----------------------------------------------------------
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  
   IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
   -----------------------------------------------------------
  </pre>
*/

public class RandomDouble implements Serializable{

    private static final int N =  624;
    private static final int M = 397;
    private static final long MATRIX_A = 0x9908b0dfL; // constant vector a
    private static final long UMASK = 0x80000000L; 
    // most significant w-r bits
    private static final long LMASK = 0x7fffffffL; 
    // least significant r bits 
    private long[] state; // the array for the state vector
    private int left = 1;
    private int initf = 0;
    private int selectedIndex;



    /** Constructor to initialize state[N] with a seed */
    public RandomDouble(long s) {
	int j;
	state =  new long[N];
	state[0]= s & 0xffffffffL;

	for (j=1; j<N; j++) {
	    state[j] = (1812433253L * (state[j-1] ^ (state[j-1] >> 30)) + (long )j); 
	    // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
	    // In the previous versions, MSBs of the seed affect   
	    // only MSBs of the array state[].                     
	    // 2002/01/09 modified by Makoto Matsumoto             
	    state[j] &= 0xffffffffL;  // for >32 bit machines 
	}
	left = 1; initf = 1;
    }

    public long MixBits(long u, long v){
	return ( ((u) & UMASK) | ((v) & LMASK) );
    }

    public long Twist(long u, long v){
	return ((MixBits(u,v) >> 1) ^ (((v) &1L)>0L ? MATRIX_A : 0L));
    }


    public void nextState(){
	int j;
	int index = 0;

	left = N;
	selectedIndex = 0;
    
	for (j=N-M+1; (--j)>0; index++) 
	    state[index] = 
		state[M+index] ^ Twist(state[index], state[index+1]);

	for (j=M; (--j)>0; index++) 
	    state[index] = 
		state[M-N+index] ^ Twist(state[index], state[index+1]);

	state[index] = state[M-N+index] ^ Twist(state[index], state[0]);
    }

    /** generates a random number on (0,1)-real-interval */
    public double nextDouble(){
	long y;

	if (--left == 0) nextState();
	y = state[selectedIndex];

	// Tempering
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680L;
	y ^= (y << 15) & 0xefc60000L;
	y ^= (y >> 18);

	selectedIndex++;

	return (((double )y + 0.5) * (1.0/4294967296.0)); 
	// divided by 2^32 
    }

}




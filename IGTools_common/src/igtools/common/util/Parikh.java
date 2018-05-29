package igtools.common.util;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;

public class Parikh {

	
	public static void getParikh(String line, Vector4Int vect){
		for(int i=0; i<vect.values.length; i++){
			vect.values[i] = 0;
		}
		B3LLSequence b3ll = new B3LLSequence(line);
		for(int i=0; i<b3ll.length(); i++)
			vect.values[b3ll.getB3(i)]++;
	}
	
	public static void getParikh(B3Nucleotide[] word, Vector4Int vect){
		for(int i=0; i<vect.values.length; i++){
			vect.values[i] = 0;
		}

		for(int i=0; i<word.length; i++)
			vect.values[word[i].code()]++;
	}
	
	
	public static String toString(Vector4Int vect){
		return (
				"A["+vect.values[0]+"]" + 
				"C["+vect.values[1]+"]" +
				"G["+vect.values[2]+"]" + 
				"T["+vect.values[3]+"]");
	}
	
	public static long nofAnagrams(int[] values, int k){
		long den = 1;
		for(int i=0; i<values.length; i++)
			den *= Maths.factorial(values[i]);
		return Maths.factorial(k) / den;
	}
        
    public static long nof4Anagrams(int[] values){
		long den = 1;
		for(int i=0; i<values.length; i++)
			den *= Maths.factorial(values[i]);
		return Maths.factorial( values[0] + values[1] + values[2] + values[3] ) / den;
	}
    
    public static long nofParikhClasses(int k){
    	return (Maths.factorial(3+k))/(Maths.factorial(k) * Maths.factorial(3));
    }
}

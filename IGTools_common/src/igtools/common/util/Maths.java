package igtools.common.util;

public class Maths {
	
//	public static long factorial(long k){
//		if(k == 0) return 1;
//		return k * factorial(k-1);
//	}
	
	public static long factorial(long k){
		long f = 1;
		for(long i=2; i<=k; i++)
			f *= i;
		return f;
	}
	
	public static long nofAnagrams(int[] values, int k){
		long den = 1;
		for(int i=0; i<values.length; i++)
			den *= factorial(values[i]);
		return factorial(k) / den;
	}
        
        public static long nof4Anagrams(int[] values){
		long den = 1;
		for(int i=0; i<values.length; i++)
			den *= factorial(values[i]);
		return factorial( values[0] + values[1] + values[2] + values[3] ) / den;
	}
}

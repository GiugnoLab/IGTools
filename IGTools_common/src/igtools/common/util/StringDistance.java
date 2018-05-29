package igtools.common.util;

public class StringDistance {

	public static enum DistanceType{
		Levenshtein,
		DamerauLevenshtein,
		Hamming,
		JaroWinkler,
		SmithWaterman,
		NGram,
		MarkovChain,
		CommonPrefixLength
	}
	
	public static int distance(String a, String b, DistanceType dtype){		
		switch(dtype){
		case Levenshtein:
			return levenshtein(a, b);
		case DamerauLevenshtein:
			return damerauLevenshtein(a, b);
		case Hamming:
			return hamming(a, b);
		case JaroWinkler:
			return jaroWinkler(a,b);
		case SmithWaterman:
			break;
		case NGram:
			break;
		case MarkovChain:
			break;
		case CommonPrefixLength:
			return commonPrefixLength(a, b);
		}
		
		return -1;
	}
	
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static int levenshtein(String a, String b){
		
		int[][] matrix = new int[a.length()+1][b.length()+1];

		for (int i = 0; i <= a.length(); i++){
			matrix[i][0] = i;
		}
		for (int j = 0; j <= b.length(); j++){
			matrix[0][j] = j;
		}
		
		for (int i = 1; i < matrix.length; i++){
			for (int j = 1; j < matrix[i].length; j++){
				if (a.charAt(i-1) == b.charAt(j-1)){
					matrix[i][j] = matrix[i-1][j-1];
				}
				else{
					int minimum = Integer.MAX_VALUE;
					if ((matrix[i-1][j])+1 < minimum){
						minimum = (matrix[i-1][j])+1;
					}

					if ((matrix[i][j-1])+1 < minimum){
						minimum = (matrix[i][j-1])+1;
					}

					if ((matrix[i-1][j-1])+1 < minimum){
						minimum = (matrix[i-1][j-1])+1;
					}

					matrix[i][j] = minimum;
				}
			}
		}
		
		return matrix[a.length()][b.length()];
	}
	
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static int damerauLevenshtein(String a, String b){
		if(a.length()>0 && b.length()>0){
			
			int[][] matrix = new int[a.length()+1][b.length()+1];
			
			int cost = -1;
			int del, sub, ins;

			matrix = new int[a.length()+1][b.length()+1];

			for (int i = 0; i <= a.length(); i++){
				matrix[i][0] = i;
			}
			for (int i = 0; i <= b.length(); i++){
				matrix[0][i] = i;
			}

			for (int i = 1; i <= a.length(); i++){
				for (int j = 1; j <= b.length(); j++){
					if (a.charAt(i-1) == b.charAt(j-1)){
						cost = 0;
					}
					else{
						cost = 1;
					}

					del = matrix[i-1][j]+1;
					ins = matrix[i][j-1]+1;
					sub = matrix[i-1][j-1]+cost;

					matrix[i][j] = damerauLevenshtein_minimum(del,ins,sub);

					if ((i > 1) && (j > 1) && (a.charAt(i-1) == b.charAt(j-2)) && (a.charAt(i-2) == b.charAt(j-1))){
						matrix[i][j] = damerauLevenshtein_minimum(matrix[i][j], matrix[i-2][j-2]+cost);
					}
				}
			}
			
			return matrix[a.length()+1][b.length()+1];
		}
		return -1;
	}
	private static int damerauLevenshtein_minimum(int d, int i, int s){
		int m = Integer.MAX_VALUE;
		if (d < m) m = d;
		if (i < m) m = i;
		if (s < m) m = s;
		return m;
	}

	private static int damerauLevenshtein_minimum(int d, int t){
		int m = Integer.MAX_VALUE;
		if (d < m) m = d;
		if (t < m) m = t;

		return m;
	}
	
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static int hamming(String a, String b){
		int l = (a.length() <= b.length()) ? a.length() : b.length();
		int m = (a.length() >= b.length()) ? a.length() : b.length();
		int dist = 0;
		for(int i=0; i<l; i++){
			if(a.charAt(i) != b.charAt(i)) dist++;
		}
		dist += m-l;
		return dist;
	}
	
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static int jaroWinkler(String a, String b){
		String theMatchA = "";
		String theMatchB = "";
		int mRange = -1;
		
		mRange = Math.max(a.length(), b.length()) / 2 - 1;

		double res = -1;

		int m = jaroWinkler_getMatch(a,b,mRange, theMatchA, theMatchB);
		int t = 0;
		if (jaroWinkler_getMissMatch(mRange, theMatchA, theMatchB) > 0){
			t = (jaroWinkler_getMissMatch(mRange, theMatchA, theMatchB) / jaroWinkler_getMissMatch(mRange, theMatchA, theMatchB));
		}

		int l1 = a.length();
		int l2 = b.length();

		double f = 0.3333;
		double mt = (double)(m-t)/m;
		double jw = f * ((double)m/l1+(double)m/l2+(double)mt);
		res = jw + jaroWinkler_getCommonPrefix(a,b) * (0.1*(1.0 - jw));

		return (int) res;
		
	}
	
	private static int jaroWinkler_getMatch(String a, String b, int mRange, String theMatchA, String theMatchB){
		theMatchA = "";
		theMatchB = "";

		int matches = 0;

		for (int i = 0; i < a.length(); i++){
			//Look backward
			int counter = 0;
			while(counter <= mRange && i >= 0 && counter <= i){
				if (a.charAt(i) == b.charAt(i - counter)){
					matches++;
					theMatchA = theMatchA + a.charAt(i);
					theMatchB = theMatchB + b.charAt(i);
				}
				counter++;
			}

			//Look forward
			counter = 1;
			while(counter <= mRange && i < b.length() && counter + i < b.length()){
				if (a.charAt(i) == b.charAt(i + counter)){
					matches++;
					theMatchA = theMatchA + a.charAt(i);
					theMatchB = theMatchB + b.charAt(i);
				}
				counter++;
			}
		}
		return matches;
	}

	private static int jaroWinkler_getMissMatch(int mRange, String theMatchA, String theMatchB){
		int transPositions = 0;

		for (int i = 0; i < theMatchA.length(); i++){
			//Look Backward
			int counter = 0;
			while(counter <= mRange && i >= 0 && counter <= i){
				if (theMatchA.charAt(i) == theMatchB.charAt(i - counter) && counter > 0){
					transPositions++;
				}
				counter++;
			}

			//Look forward
			counter = 1;
			while(counter <= mRange && i < theMatchB.length() && (counter + i) < theMatchB.length()){
				if (theMatchA.charAt(i) == theMatchB.charAt(i + counter) && counter > 0){
					transPositions++;
				}
				counter++;
			}
		}
		return transPositions;
	}

	private static int jaroWinkler_getCommonPrefix(String a, String b){
		int cp = 0;
		for (int i = 0; i < 4; i++){
			if (a.charAt(i) == b.charAt(i)) cp++;
		}
		return cp;
	}
	
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static int commonPrefixLength(String a, String b){
		int score = 0;
		for(int i=0; i<a.length() && i<b.length(); i++){
			if(a.charAt(i) == b.charAt(i))
				score++;
			else break;
		}
		return score;
	}
	
}

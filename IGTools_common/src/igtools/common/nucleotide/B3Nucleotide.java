package igtools.common.nucleotide;

import java.util.Comparator;


public enum B3Nucleotide implements INucleotide{
	A	((byte)0), 
	C	((byte)1), 
	G	((byte)2), 
	T	((byte)3), 
	N 	((byte)4),
	NULL ((byte)5);
	
	public static final  byte A_CODE = (byte)0;
	public static final byte C_CODE = (byte)1;
	public static final byte G_CODE = (byte)2;
	public static final byte T_CODE = (byte)3;
	public static final byte N_CODE = (byte)4;
	public static final byte NULL_CODE = (byte)5;
	
	private byte code;
	
	private B3Nucleotide(byte c){
		code = c;
	}
	
	@Override
	public byte code(){
		return code;
	}
	
	public char getChar(){
		switch(code){
		case A_CODE:
			return 'A';
		case C_CODE:
			return 'C';
		case G_CODE:
			return 'G';
		case T_CODE:
			return 'T';
		case N_CODE:
			return 'N';
		}
		//return B3Nucleotide.NULL;
		return '#';
	}
	
	public static B3Nucleotide by(int i){
		switch(i){
		case 0:
			return B3Nucleotide.A;
		case 1:
			return B3Nucleotide.C;
		case 2:
			return B3Nucleotide.G;
		case 3:
			return B3Nucleotide.T;
		case 4:
			return B3Nucleotide.N;
		}
		//return B3Nucleotide.NULL;
		return null;
	}
	
	public static B3Nucleotide by(char c){
		switch(c){
		case 'A':
		case 'a':
			return B3Nucleotide.A;
		case 'C':
		case 'c':
			return B3Nucleotide.C;
		case 'G':
		case 'g':
			return B3Nucleotide.G;
		case 'T':
		case 't':
			return B3Nucleotide.T;
		case 'U':
		case 'u':
		case 'R':
		case 'r':
		case 'Y':
		case 'y':
		case 'K':
		case 'k':
		case 'M':
		case 'm':
		case 'S':
		case 's':
		case 'W':
		case 'w':
		case 'B':
		case 'b':
		case 'D':
		case 'd':
		case 'H':
		case 'h':
		case 'V':
		case 'v':
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		case '-':
			return B3Nucleotide.N;
		}
		//return B3Nucleotide.NULL;
		return null;
	}
	
	public static char charFor(int code){
		switch(code){
		case A_CODE:
			return 'A';
		case C_CODE:
			return 'C';
		case G_CODE:
			return 'G';
		case T_CODE:
			return 'T';
		case N_CODE:
			return 'N';
		}
		//return B3Nucleotide.NULL;
		return '#';
	}
	
	public static int codeFor(char c){
		switch(c){
		case 'A':
		case 'a':
			return A_CODE;
		case 'C':
		case 'c':
			return C_CODE;
		case 'G':
		case 'g':
			return G_CODE;
		case 'T':
		case 't':
			return T_CODE;
		case 'U':
		case 'u':
		case 'R':
		case 'r':
		case 'Y':
		case 'y':
		case 'K':
		case 'k':
		case 'M':
		case 'm':
		case 'S':
		case 's':
		case 'W':
		case 'w':
		case 'B':
		case 'b':
		case 'D':
		case 'd':
		case 'H':
		case 'h':
		case 'V':
		case 'v':
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		case '-':
			return N_CODE;
		default:
			return B3Nucleotide.NULL_CODE;
		}
		//return B3Nucleotide.NULL_CODE;
	}
	
	public static boolean isBAD(int code){
		return code == N_CODE;
	}
        
        
        public static void copy(B3Nucleotide[] source, B3Nucleotide[] target){
        	for(int i=0; i<source.length; i++){
        		target[i] = source[i];
        	}
        }
        
        
        public static String toString(B3Nucleotide[] kmer){
            String s = "";
            for(int i=0; i<kmer.length; i++)
                s += kmer[i];
            return s;
        }
        public static String toString(B3Nucleotide[] kmer, int length){
            String s = "";
            for(int i=0; i<kmer.length && i<length; i++)
                s += kmer[i];
            return s;
        }
        
        public static String toString(B3Nucleotide[] kmer, int start, int end){
            String s = "";
            for(int i=start; i<end; i++)
                s += kmer[i];
            return s;
        }
        
        
        public static Comparator<B3Nucleotide[]> getComparator(){
            return new Comparator<B3Nucleotide[]>(){
                @Override
                public int compare(B3Nucleotide[] o1, B3Nucleotide[] o2) {
                    for(int i=0; i<(o1.length < o2.length ? o1.length : o2.length ); i++){
                        if(o1[i] != o2[i]){
                            return o1[i].code() - o2[i].code();
                        }
                    }
                    return 0;
                }
            };
        }
        
        
        public static B3Nucleotide[] toB3(String s){
        	B3Nucleotide[] ns = new B3Nucleotide[s.length()];
        	for(int i=0; i<s.length(); i++){
        		ns[i] = B3Nucleotide.by(s.charAt(i));
        	}
        	return ns;
        }
        
        public static void toB3(String s, B3Nucleotide[] ns){
        	for(int i=0; i<s.length(); i++){
        		ns[i] = B3Nucleotide.by(s.charAt(i));
        	}
        }
        
        public static B3Nucleotide[] toB3(B2Nucleotide[] b2kmer){
            B3Nucleotide[] ns = new B3Nucleotide[b2kmer.length];
            
            for(int i=0;i<ns.length; i++){
                ns[i] = B3Nucleotide.by(B2Nucleotide.charFor(b2kmer[i].code()));
            }
            
            return ns;
        }
        
        public static int compare(B3Nucleotide[] n1, B3Nucleotide[] n2){
//        	System.out.println("compare "+toString(n1) +","+toString(n2));
        	for(int i=0; i<n1.length; i++){
        		if(n1[i].code != n2[i].code)
        			return n1[i].code - n2[i].code;
        	}
        	return 0;
        }
        
        public static int compareToC(B3Nucleotide[] ns){
        	for(int i=0; i<ns.length; i++){
        		if(ns[i] != rcs[ns[i].code]){
        			return ns[i].code - rcs_code[ns[i].code];
        		}
        	}
        	return 0;
        }
        
        public static int compareToRC(B3Nucleotide[] ns){
        	for(int i=0; i<ns.length; i++){
        		if(ns[i] != rcs[ns[ns.length - (i+1)].code]){
        			return ns[i].code - rcs_code[ns[ns.length - (i+1)].code];
        		}
        	}
        	return 0;
        }
        
        public static boolean areEqual(B3Nucleotide[] n1, B3Nucleotide[] n2){
//        	System.out.println("equals "+toString(n1) +","+toString(n2));
        	for(int i=0; i<n1.length; i++){
        		if(n1[i].code != n2[i].code)
        			return false;
        	}
        	return true;
        }
        
//        public static B3Nucleotide[] toB3(String kmer){
//            B3Nucleotide[] ns = new B3Nucleotide[kmer.length()];
//            
//            
//            return ns;
//        }
        
        public static int compl(int code){
        	return rcs_code[code];
        }
        public static B3Nucleotide complB3(int code){
        	return rcs[code];
        }
        public static B3Nucleotide compl(B3Nucleotide n){
        	return rcs[n.code];
        }
        public static B3Nucleotide complement(B3Nucleotide n){
        	switch(n){
			case A:
				return B3Nucleotide.T;
			case C:
				return B3Nucleotide.G;
			case G:
				return B3Nucleotide.C;
			case N:
				return B3Nucleotide.N;
			case NULL:
				return B3Nucleotide.NULL;
			case T:
				return B3Nucleotide.A;
			default:
				return B3Nucleotide.NULL;
        	}
        }
        
        
        
//        public static final  byte A_CODE = (byte)0;
//    	public static final byte C_CODE = (byte)1;
//    	public static final byte G_CODE = (byte)2;
//    	public static final byte T_CODE = (byte)3;
//    	public static final byte N_CODE = (byte)4;
//    	public static final byte NULL_CODE = (byte)5;
        public static final B3Nucleotide[] rcs = {B3Nucleotide.T, B3Nucleotide.G, B3Nucleotide.C, B3Nucleotide.A, B3Nucleotide.N, B3Nucleotide.NULL};
        public static final int[] rcs_code = {B3Nucleotide.T.code, B3Nucleotide.G.code, B3Nucleotide.C.code, B3Nucleotide.A.code, B3Nucleotide.N.code, B3Nucleotide.NULL.code};
        
        public static void complement(B3Nucleotide[] ns){
        	for(int i=0; i<ns.length; i++){
        		//ns[i] = rc(ns[i]);
        		ns[i] = rcs[ns[i].code];
        	}
        }
        
        public static void reverseComplement(B3Nucleotide[] ns){
        	//final int to = ((int)Math.ceil(ns.length / 2)) + 1;
        	int to = (ns.length / 2);
        	if(ns.length%2 != 0)
        		to++;
        	
        	//System.out.println(ns.length + " "+to);
        	B3Nucleotide tmp;
        	for(int i=0; i<to; i++){
        		
        		//System.out.println(i+" -> "+(ns.length - i - 1));
        		
        		tmp = rcs[ns[i].code];
        		ns[i] = rcs[ ns[ns.length - i - 1].code ];
        		ns[ns.length - i - 1] = tmp;
        	}
        }
        public static B3Nucleotide[] reverseComplementCopy(B3Nucleotide[] ns){
        	B3Nucleotide[] rc = new B3Nucleotide[ns.length];
        	copy(ns, rc);
        	reverseComplement(rc);
        	return rc;
        }
        
        public static void reverseComplement(B3Nucleotide[] ns, B3Nucleotide[] rc){
        	copy(ns, rc);
        	reverseComplement(rc);
        }
        
        public static long toLexicoOrder(B3Nucleotide[] ns){
        	long ret = 0;
        	
        	long p4 = 1;
        	for(int i=0; i<ns.length - 1; i++){
        		p4 *= 4;
        	}
        	
        	for(int i=0; i<ns.length; i++){
        		ret += ns[i].code *  p4;
        		p4 = p4 / 4;
        	}
        	
        	return ret;
        }
        
        public static long toLexicoOrder(B3Nucleotide[] ns, int start, int end){
        	long ret = 0;
        	
        	long p4 = 1;
        	for(int i=start; i<end - 1; i++){
        		p4 *= 4;
        	}
        	
        	for(int i=start; i<end; i++){
        		ret += ns[i].code *  p4;
        		p4 = p4 / 4;
        	}
        	
        	return ret;
        }
        
}

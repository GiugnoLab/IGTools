package igtools.common.nucleotide;



public enum B2Nucleotide implements INucleotide{
	
	A	((byte)0), 
	C	((byte)1), 
	G	((byte)2), 
	T	((byte)3);
        
    public static final  byte A_CODE = (byte)0;
	public static final byte C_CODE = (byte)1;
	public static final byte G_CODE = (byte)2;
	public static final byte T_CODE = (byte)3;
        

	private byte code;

	private B2Nucleotide(byte c){
			code = c;
	}

	@Override
	public byte code(){
			return code;
	}

	public static B2Nucleotide by(int i){
		switch(i){
		case 0:
		return B2Nucleotide.A;
		case 1:
		return B2Nucleotide.C;
		case 2:
		return B2Nucleotide.G;
		case 3:
		return B2Nucleotide.T;
		}
		return null;
	}

	public static B2Nucleotide by(char c){
		switch(c){
			case 'A':
			case 'a':
			return B2Nucleotide.A;
			case 'C':
			case 'c':
			return B2Nucleotide.C;
			case 'G':
			case 'g':
			return B2Nucleotide.G;
			case 'T':
			case 't':
			return B2Nucleotide.T;
		}
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
		}
		//return B3Nucleotide.NULL;
		return '#';
	}
        
        
        public static String toString(B2Nucleotide[] kmer){
            String s = "";
            for(int i=0; i<kmer.length; i++)
                s += kmer[i];
            return s;
        }
}

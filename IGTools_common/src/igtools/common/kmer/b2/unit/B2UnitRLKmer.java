package igtools.common.kmer.b2.unit;

import igtools.common.nucleotide.B2Nucleotide;




/**
 * Right Most, Left to Right, Integer, Kmer (max k = 16 for 2bits representation)
 * 
 * [00...00|kmer]
 * 
 * 
 * @author bovi
 *
 */
public class B2UnitRLKmer {
	public static final int BITS_PER_MER = 2;
	
	public static final int  MAX_K = 16;
	
	public static final int[] MAX_CODE;
	public static final int[] MAX_SHIFT;
	
	
	static{
		MAX_CODE = new int[MAX_K + 1];//[0...16]
		MAX_CODE[0] = 0x0;
		int mask = 0x0;
		int shift = 0;
		for(int i=1; i<17; i++){
			mask |= 0x3 << shift;
			MAX_CODE[i] = mask;
			shift += 2;
		}
		
		MAX_SHIFT = new int[MAX_K + 1];//[0...16]
		MAX_SHIFT[0] = 0x0;
		shift = 0;
		for(int i=1; i<17; i++){
			MAX_SHIFT[i] = shift;
			shift += 2;
		}
	}
	
	
	
	private int code = 0x0;
	
	public B2UnitRLKmer(int code){
		this.code = code;
	}
	
	
	/**
	 * 
	 * @param ns must be in Left to Right form
	 */
	public B2UnitRLKmer(B2Nucleotide[] ns) throws Exception{
		if(ns.length > MAX_K ) throw new Exception("Too long nucleotide sequence");
//		else if(ns.length < 1) code = 0x0;
//		else{
			int shift = ((ns.length - 1) * 2);
			for(int i=0; i<ns.length; i++){
				code |= ns[i].code() << shift;
				shift -= 2;
			}
//		}
	}
	
	public int code(){
		return code;
	}
        public void setCode(int code){
            this.code = code;
        }
	
	/**
	 * Left to Right
	 * 
	 * @return
	 */
	public B2Nucleotide[] nucleotides(){
		return nucleotides(16);
	}
	
	/**
	 * Left to Right
	 * 
	 * @param size
	 * @return
	 */
	public B2Nucleotide[] nucleotides(int size){
		if(size > MAX_K) size = MAX_K;
		if(size <1) size = 1;
		B2Nucleotide[] ns = new B2Nucleotide[size];
		int shift = ((size - 1) * 2);//30 for k = 16
		for(int i=0; i<size; i++){
			ns[i] =  B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return ns;
	}
	
	
	/**
	 * Left to Right
	 */
	public String toString(){
		String s = "";
		int shift = ((MAX_K - 1) * 2);
		for(int i=0; i<MAX_K; i++){
			s += B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return s;
	}
        
        public String toString(int k){
		String s = "";
		int shift = ((k - 1) * 2);
		for(int i=0; i<k; i++){
			s += B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return s;
	}

}

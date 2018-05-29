package igtools.common.kmer.b2.unit;

import igtools.common.nucleotide.B2Nucleotide;



/**
 * -> Left to Right  12mer 
 * 
 *  First nucleotide in position 22, second one in position 20... and so on
 *  
 *  32bit word
 *  
 *  nucl	 ____4___ ________12_________
 *  bits 	 ____8___ ________24_________
 *  		|00000000|First..........Last|
 *  		|not used| 12mer ->	L to R	 |
 * 
 * @author bovi
 *
 */
public class B2UnitLR12mer {
	public static final int MAX_VALUE = 0xFFFFFF;//1 * 24bit;  12 * G 
	public static final int FIRST_SHIFT = 22;//bits
	public static final int LENGTH = 12;//nucleotides
	
	private int code = 0x0;
	
	public B2UnitLR12mer(int kmer){
		this.code = kmer & MAX_VALUE;
	}
	
	
	/**
	 * 
	 * @param ns must be in -> Left to Right format
	 * 
	 * @throws Exception is length of ns != 12
	 */
	public B2UnitLR12mer(B2Nucleotide[] ns) throws Exception{
		if(ns.length != 12) throw new Exception("Too long nucleotide sequence");
		
		int shift = FIRST_SHIFT;
		for(int i=0; i<ns.length; i++){
			code |= ns[i].code() << shift;
			shift -= 2;
		}
	}

	public int code(){
		return code;
	}
	
	
	public B2Nucleotide[] nucleotides(){
		B2Nucleotide[] ns = new B2Nucleotide[LENGTH];
		int shift = FIRST_SHIFT;
		for(int i=0; i<LENGTH; i++){
			ns[i] =  B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return ns;
	}
	
	public B2Nucleotide[] nucleotides(int size){
		if(size > LENGTH) size = LENGTH;
		if(size <0) size = 0;
		B2Nucleotide[] ns = new B2Nucleotide[size];
		int shift = FIRST_SHIFT;
		for(int i=0; i<size; i++){
			ns[i] =  B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return ns;
	}
	
	public String toString(){
		String s = "";
		int shift = FIRST_SHIFT;
		for(int i=0; i<LENGTH; i++){
			s += B2Nucleotide.by((code >>> shift) & 0x3);
			shift -= 2;
		}
		return s;
	}
}

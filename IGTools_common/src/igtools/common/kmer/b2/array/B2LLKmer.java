package igtools.common.kmer.b2.array;

import igtools.common.kmer.IArrayKmer;
import igtools.common.nucleotide.B2Nucleotide;


/**
 * Left Most, Left to Right, generic kmer.
 * 
 * First nucleotide in position 30 of first array element
 * 
 *  0      ... n
 *  [km...]...[...er|00...00]
 * 
 * @author bovi
 *
 */
public class B2LLKmer implements IArrayKmer, Comparable<B2LLKmer>{
	public final static int BITS_PER_MER = 2;
	public final static int FIRST_SHIFT = 30;
	public final static int BITS_PER_BLOCK = 32;
	public final static int KMERS_PER_BLOCK = BITS_PER_BLOCK / BITS_PER_MER;
	
	private int length;
	private int[] kmer = null;

	
	public static int dataBinsFor(int k){
		return (int)Math.ceil((double)k / KMERS_PER_BLOCK);
	}
	
	
	public B2LLKmer(int length, int[] kmer){
		this.length = length;
		this.kmer = kmer;
	}
	public B2LLKmer(B2Nucleotide[] ns) {
		this.length = ns.length;
		kmer = new int[ (int)Math.ceil((double)ns.length / KMERS_PER_BLOCK) ];
		getFrom(ns);
	}
	
	public B2LLKmer(String s){
		this.length = s.length();
		kmer = new int[ (int)Math.ceil((double)length / KMERS_PER_BLOCK) ];
		B2Nucleotide[] ns = new B2Nucleotide[length];
		for(int i=0; i<length; i++){
			ns[i] = B2Nucleotide.by(s.charAt(i));
		}
		getFrom(ns);
	}
	
	public int length(){
		return length;
	}
	
	public int[] kmer(){
		return kmer;
	}
	
	public int[] code(){
		return kmer;
	}
	
	private void getFrom(B2Nucleotide[] ns){
		int block = 0;
		int pos = FIRST_SHIFT;
		int i = 0;
		int n;
		
		for(; i<ns.length; i++){
			if(i!=0 && ((i*BITS_PER_MER) % BITS_PER_BLOCK == 0)){
				pos = FIRST_SHIFT;
				block++;
			}
			n = ns[i].code();
			n = (n << pos);
			kmer[block] |= n;
			pos -= BITS_PER_MER;
		}
	}

	public B2Nucleotide[] nucleotides() {
		B2Nucleotide[] ns = new B2Nucleotide[this.length];

		int pos = 0;
		int mask = 0xC0000000;
		int shift = FIRST_SHIFT;
		
		for(int i=0; i<ns.length; i++){
			if(i!=0 && ((i*BITS_PER_MER) % BITS_PER_BLOCK == 0)){
				mask = 0xC0000000;
				
				pos++;
				shift = FIRST_SHIFT;
			}
			ns[i] = B2Nucleotide.by( (kmer[pos] & mask) >>> shift);
			mask = mask >>> BITS_PER_MER;
			shift -= BITS_PER_MER;
		}
		
		
		return ns;
	}
	

	@Override
	public boolean equals(Object o) {
		if(o instanceof B2LLKmer){
			if(this.length == ((B2LLKmer)o).length()){
				int[] b = ((B2LLKmer)o).kmer();
				for(int i=0; i<kmer.length; i++){
					if(kmer[i] != b[i])
						return false;
				}
				return true;
			}
		}
		return false;
	}
	
	@Override
	public int hashCode() {
		int c = 0;
		for(int i=0; i<kmer.length; i++){
			c += kmer[i];
		}
		return c;
	}
	
	
	@Override
	public String toString(){
		String s = "";
		B2Nucleotide[] ns = this.nucleotides();
		for(int i=0; i<ns.length; i++){
			s = s + ns[i];
		}
		return s;
	}
	@Override
	public int compareTo(B2LLKmer o) {
//		int ret = 0;
//		int min = kmer.length <= o.kmer.length ? kmer.length : o.kmer.length; 
//		for(int i=0; i<min; i++){
//			ret += kmer[i] - o.kmer[i];
//		}
//		return ret;
		
		
		//ODIO JAVA
		
		
		int ret = 0;
		int min = kmer.length <= o.kmer.length ? kmer.length : o.kmer.length; 
		for(int i=0; i<min; i++){
			if(kmer[i] != o.kmer[i]){
				if( ((long)(kmer[i] & 0xffffffffL)) < ((long)(o.kmer[i] & 0xffffffffL)) )
					return -1;
				else
					return 1;
			}
		}
		if(ret == 0 && kmer.length != o.kmer.length){
			if(kmer.length < o.kmer.length)
				return -1;
			else
				return 1;
		}
		
		return ret;
	}
}

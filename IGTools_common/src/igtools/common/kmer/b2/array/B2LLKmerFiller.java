package igtools.common.kmer.b2.array;

public class B2LLKmerFiller {

	private int[] code;
	private int iblock = 0x0;
	private int shift = B2LLKmer.FIRST_SHIFT;
	
	
	public B2LLKmerFiller(int length){
		code = new int[(length / B2LLKmer.BITS_PER_MER) + ( (length % B2LLKmer.BITS_PER_MER) > 0 ? 1 : 0)];
		//code = new int[ (int)Math.ceil((double)ns.length / 16) ];
	}
	
	public void push(int c){		
		code[iblock] |= c << shift;
		shift -= B2LLKmer.BITS_PER_MER;
		if(shift < 0){
			shift = B2LLKmer.FIRST_SHIFT;
			iblock++;
		}
	}
	
	
	public int[] code(){
		return code;
	}
}

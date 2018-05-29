package igtools.common.kmer.b3.unit;

public class B3UnitLLKmer {
	
	public final static int BITS_PER_BLOCK = 32;
	public final static int BITS_PER_MER = 3;
	public final static int MERS_PER_BLOCK = BITS_PER_BLOCK / BITS_PER_MER;
	public final static int FIRST_SHIFT = BITS_PER_BLOCK - BITS_PER_MER;//32 - 3 = 29 
	public final static int MASK = 0x7;//111 
	public final static int ALL_ONES = 0xFFFFFFFF;
	

}

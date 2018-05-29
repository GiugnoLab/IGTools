package igtools.common.sequence;

import igtools.common.kmer.b2.array.B2LLKmer;
import igtools.common.nucleotide.B3Nucleotide;

public interface B3Sequence extends ISequence<B2LLKmer>, CharSequence{

	public void clear(int length);
	public int length();
	
	public B2LLKmer get(int pos, int length);
	public int getB2(int pos);
	
	public int getB3(int pos);
	public int getB3(int pos,B3Nucleotide[] ns);
	
	public void setB3(int pos, int code);

	public int indexOf(int code, int from, int to);
	public int lcp(int pos_a, int pos_b);
	
	public byte[] toB3ByteArray();
	public byte[] toB2ByteArray();
	
	public int countBads();
}

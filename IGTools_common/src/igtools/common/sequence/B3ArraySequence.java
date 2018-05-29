package igtools.common.sequence;

import igtools.common.kmer.b2.array.B2LLKmer;
import igtools.common.nucleotide.B3Nucleotide;

public class B3ArraySequence implements B3Sequence{
	
	private B3Nucleotide[] iseq;
	private int iseq_start;
	private int iseq_end;
	private int rlength;
	
	public B3ArraySequence(B3Nucleotide[] iseq){
		this.iseq = iseq;
		iseq_start = 0;
		iseq_end = iseq.length;
		rlength = iseq.length;
	}
	
	public B3ArraySequence(B3Nucleotide[] iseq, int start, int end){
		this.iseq = iseq;
		iseq_start = start;
		iseq_end = end;
		rlength = end-start; 
	}

	@Override
	public char charAt(int index) {
		if(index>-1 && index<rlength)
			return iseq[index - iseq_start].getChar();
		return B3Nucleotide.NULL.getChar();
	}

	@Override
	public CharSequence subSequence(int start, int end) {
		return new B3ArraySequence(iseq, iseq_start + start, iseq_start + end);
	}

	@Override
	public void clear(int length) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int length() {
		return rlength;
	}

	@Override
	public B2LLKmer get(int pos, int length) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int getB2(int pos) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int getB3(int pos) {
		return iseq[iseq_start + pos].code();
	}

	@Override
	public int getB3(int pos, B3Nucleotide[] ns) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public void setB3(int pos, int code) {
		iseq[iseq_start + pos] = B3Nucleotide.by(code);
	}

	@Override
	public int indexOf(int code, int from, int to) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int lcp(int pos_a, int pos_b) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public byte[] toB3ByteArray() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public byte[] toB2ByteArray() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int countBads() {
		int res = 0;
		for(int i=iseq_start; i<iseq_end; i++){
			if(iseq[i] == B3Nucleotide.N)
				res++;
		}
		return res;
	}

}

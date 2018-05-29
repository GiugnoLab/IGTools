package igtools.common.sequence;

import igtools.common.kmer.b2.array.B2LLKmer;
import igtools.common.nucleotide.B3Nucleotide;


//TODO
public class OnDiskB3LLSequence implements B3Sequence {
	
	
	
	
	
	public OnDiskB3LLSequence(String file){
		
	}
	

	@Override
	public char charAt(int index) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public CharSequence subSequence(int start, int end) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void clear(int length) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int length() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public B2LLKmer get(int pos, int length) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getB2(int pos) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getB3(int pos) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getB3(int pos, B3Nucleotide[] ns) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setB3(int pos, int code) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int indexOf(int code, int from, int to) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public byte[] toB3ByteArray() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte[] toB2ByteArray() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public int lcp(int pos_a, int pos_b) {
		// TODO Auto-generated method stub
		return 0;
	}


	@Override
	public int countBads() {
		// TODO Auto-generated method stub
		return 0;
	}

	
}

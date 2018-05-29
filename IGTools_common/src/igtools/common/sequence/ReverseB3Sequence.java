package igtools.common.sequence;

import igtools.common.kmer.b2.array.B2LLKmer;
import igtools.common.nucleotide.B3Nucleotide;

public class ReverseB3Sequence implements B3Sequence{

	
	private B3Sequence base;
	private int l,r;
	
	public ReverseB3Sequence(B3Sequence _base, int start, int end){
		base = _base;
		l = start;
		r = end;
	}
	public ReverseB3Sequence(B3Sequence _base){
		this(_base, 0, _base.length());
	}
	
	
	@Override
	public char charAt(int index) {
		return base.charAt(r - index -1);
	}
	@Override
	public B3Sequence subSequence(int start, int end) {
		return new ReverseB3Sequence(base, r-end, r-start);
	}
	@Override
	public void clear(int length) {
		base.clear(length);
	}
	@Override
	public int length() {
		return r-l;
	}
	@Override
	public B2LLKmer get(int pos, int length) {
		// TODO Auto-generated method stub
		throw new Error("Unstatiated method");
		//return null;
	}
	@Override
	public int getB2(int pos) {
		// TODO Auto-generated method stub
//		return 0;
		throw new Error("Unstatiated method");
	}
	@Override
	public int getB3(int pos) {
		return base.getB3(r -pos -1);
	}
	@Override
	public int getB3(int pos, B3Nucleotide[] ns) {
		// TODO Auto-generated method stub
//		return 0;
		int i = 0;
		for(; i<ns.length && pos+i< this.length(); i++)
			ns[i] = B3Nucleotide.by(getB3(pos + i));
		return i;
	}
	@Override
	public void setB3(int pos, int code) {
		// TODO Auto-generated method stub
		//base.setB3(r -pos -1, code);
		throw new Error("Unstatiated method");
	}
	
	public String toString(){
		//String s  = "s:"+this.length()+":";
		String s  = "";
		for(int i=0; i<this.length(); i++){
			s = s + this.charAt(i);
		}
		return s;
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
	//@Override
	public int lcp(int pos_a, int pos_b) {
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public int countBads() {
		// TODO Auto-generated method stub
		throw new Error("Unstatiated method");
		//return 0;
	}
}

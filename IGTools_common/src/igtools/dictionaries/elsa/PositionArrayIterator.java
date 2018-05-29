package igtools.dictionaries.elsa;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;


/**
 * Iterate over k-mers defined by their start positions.
 * IT does not check for duplicates. 
 * 
 * @author vbonnici
 *
 */
public class PositionArrayIterator implements IELSAIterator {

	
	
	private int[] p_pos;
	private int p_i;
	
	private int k;
	private B3Sequence seq;
	
	
	public PositionArrayIterator(int[] p_pos, int k, B3Sequence seq){
		this.p_pos = p_pos;
		this.k = k;
		this.seq = seq;
		this.p_i = -1;
	}
	
	
	@Override
	public IELSA elsa() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public int istart() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public int iend() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public int k() {
		return this.k;
	}
	@Override
	public B3Nucleotide[] kmer() {
		B3Nucleotide[] kmer = new B3Nucleotide[k];
		seq.getB3(p_pos[p_i], kmer);
		return kmer;
	}
	@Override
	public void kmer(B3Nucleotide[] ns) {
		seq.getB3(p_pos[p_i], ns);
	}
	@Override
	public int multiplicity() {
		return 1;
	}
	@Override
	public int[] positions() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public int[] sortedPositions() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public boolean isMinimalHapax() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public boolean isGlobalMaximalRepeat() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	@Override
	public boolean next() {
		if(this.p_i  <  this.p_pos.length - 1){
			this.p_i++;
			return true;
		}
		return false;
	}
	@Override
	public boolean prev() {
		if(this.p_i > 0){
			this.p_i--;
			return true;
		}
		return false;
	}
	@Override
	public boolean good() {
		return (this.p_i > -1) && (this.p_i < this.p_pos.length);
	}
	@Override
	public boolean hasNext() {
		return (this.p_i  <  this.p_pos.length - 1);
	}
	@Override
	public boolean hasPrev() {
		return this.p_i > 0;
	}
	@Override
	public int compare(IELSAIterator it) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	
	
	@Override
	public IELSAIterator clone(){
		return new PositionArrayIterator(this.p_pos, this.k, this.seq);
	}
	
	
}

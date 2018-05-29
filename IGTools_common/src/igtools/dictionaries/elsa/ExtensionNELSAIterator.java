package igtools.dictionaries.elsa;

import java.util.Arrays;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;

public class ExtensionNELSAIterator implements IELSAIterator{
	
	private IELSAIterator parent;
	private NELSA nelsa;
	private B3Sequence s;
	private int[] sa;
	private int[] lcp;
	private int[] ns;
	private int k;
	private int istart;//inclusive
	private int iend;//exclusive
	private int ilimit;
	private int limit;
	
	
	public ExtensionNELSAIterator(NELSA nelsa, IELSAIterator it){
		parent = it;
		this.nelsa = nelsa;
		s = nelsa.b3seq();
		
		sa = nelsa.sa();
		lcp = nelsa.lcp();
		ns = nelsa.ns();
		
		k = it.k() + 1;
		
		istart = it.istart()-1;
		iend = it.istart();
		ilimit = it.iend();
		
		limit = nelsa.length() - k;
	}
	
	
	

	@Override
	public IELSA elsa() {
		return nelsa;
	}

	@Override
	public int istart() {
		return istart;
	}

	@Override
	public int iend() {
		return iend;
	}

	@Override
	public int k() {
		return k;
	}

	@Override
	public B3Nucleotide[] kmer() {
		B3Nucleotide[] seq = new B3Nucleotide[k];
		s.getB3(sa[istart], seq);
		return seq;
	}

	@Override
	public void kmer(B3Nucleotide[] ns) {
		s.getB3(sa[istart], ns);
	}
	
	
	public int lastNucleotide(){
		return s.getB3(sa[istart] + k - 1);
	}

	@Override
	public int multiplicity() {
		return iend - istart;
	}

	@Override
	public int[] positions() {
		return Arrays.copyOfRange(sa, istart, iend);
	}

	@Override
	public int[] sortedPositions() {
		int[] pos =  Arrays.copyOfRange(sa, istart, iend);
		Arrays.sort(pos);
		return pos;
	}

	@Override
	public boolean isMinimalHapax() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean isGlobalMaximalRepeat() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean next() {
		int l = iend;
		//while(l<ilimit && (ns[l]<=k || sa[l]>limit)){
		while(l<ilimit && (ns[l]<k || sa[l]>limit)){
			l++;
		}
		
		if(l >= ilimit){
			istart = -1;
			iend = 0;
			return false;
		}
		
		istart = l;
		l++;
		//while(l<ilimit && lcp[l]>=k && sa[l]<limit && ns[l]>k){
		while(l<ilimit && lcp[l]>=k && sa[l]<limit && ns[l]>=k){
			l++;
		}
		iend = l;
		return true;
	}

	@Override
	public boolean prev() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean good() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean hasNext() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean hasPrev() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public int compare(IELSAIterator it) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}
	
	 @Override
     public IELSAIterator clone() {
		 ExtensionNELSAIterator c = new ExtensionNELSAIterator(nelsa, parent);
		 c.istart = istart;
		 c.iend = iend;
		 return c;
     }

}

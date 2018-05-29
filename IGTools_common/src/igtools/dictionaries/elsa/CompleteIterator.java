package igtools.dictionaries.elsa;

import igtools.common.nucleotide.B3Nucleotide;

public class CompleteIterator implements IELSAIterator{
	
	private int k;
	private int[] ikmer;
	private B3Nucleotide[] kmer;
	
	private int p = -1;
	
	public CompleteIterator(int k){
		this.k = k;
		ikmer = new int[k];
		kmer = new B3Nucleotide[k];
		for(int i=0;i<k;i++){
			ikmer[i] = 0;
			kmer[i] = B3Nucleotide.by(0);
		}
	}

	@Override
	public ELSA elsa() {
		return null;
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
		return k;
	}

	@Override
	public B3Nucleotide[] kmer() {
		return kmer;
	}

	@Override
	public void kmer(B3Nucleotide[] ns) {
		for(int i=0; i<ns.length; i++){
			ns[i] = kmer[i];
		}
	}

	@Override
	public int multiplicity() {
		throw new UnsupportedOperationException("Not supported yet.");
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
		if(p == -1){
			p++;
			return p<k;
		}
		else{
			p = k-1;
			while(p > -1){
				if(ikmer[p] == 3){
					ikmer[p] = 0;
					kmer[p] = B3Nucleotide.by(ikmer[p]);
	//				System.out.println("n: "+p+"\t"+ikmer[p]);
					p--;
				}
				else{
					ikmer[p]++;
					kmer[p] = B3Nucleotide.by(ikmer[p]);
	//				System.out.println("t: "+p+"\t"+ikmer[p]);
					return true;
				}
			}
			if(p == -1){
				p = 0;
				return false;
			}
			return true;
		}
	}

	@Override
	public boolean prev() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean good() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean hasNext() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean hasPrev() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int compare(IELSAIterator it) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	
	@Override
	public IELSAIterator clone(){
		return new CompleteIterator(k);
	}
	
	
	
	
	public static void main(String[] args){
		for(int k=1; k<18; k++){
			System.out.println(k);
			CompleteIterator it = new CompleteIterator(k);
			while(it.next()){
				System.out.println(B3Nucleotide.toString(it.kmer()));
			}
		}
	}
	

}

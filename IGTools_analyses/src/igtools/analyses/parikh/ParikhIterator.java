package igtools.analyses.parikh;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.IELSA;
import igtools.dictionaries.elsa.IELSAIterator;

public class ParikhIterator implements IELSAIterator{
	
	
	
	B3Nucleotide[] symbols;
	boolean[] choosen;
	
	int k;
	int[] its;
	int it;
	
	/**
	 * 
	 * @param values, must be a int[4], values[0] = #A, values[1] = #C, #G, #T
	 */
	public ParikhIterator(int k, int[] values){
		int sum = 0;
		for(int i=0; i<values.length; i++)
			sum += values[i];
		
		symbols = new B3Nucleotide[sum];
		choosen = new boolean[sum];
		
		for(int i=0; i<sum; i++)
			choosen[i] =false;
		
		
		int j = 0;
		for(int i=0; i<values[0]; i++, j++)
			symbols[j] = B3Nucleotide.A;
		for(int i=0; i<values[1]; i++, j++)
			symbols[j] = B3Nucleotide.C;
		for(int i=0; i<values[2]; i++, j++)
			symbols[j] = B3Nucleotide.G;
		for(int i=0; i<values[3]; i++, j++)
			symbols[j] = B3Nucleotide.T;
		
		
		this.k = k;
		its = new int[k];
		restart();
	}
	
	private void restart(){
		it = 0;
		
		for(int i=0; i<k; i++)
			its[i] = -1;
		
		for(int i=0; i<choosen.length; i++)
			choosen[i] =false;
	}
	

	@Override
	public IELSA elsa() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public int istart() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public int iend() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public int k() {
		return this.k;
	}

	@Override
	public B3Nucleotide[] kmer() {
		B3Nucleotide[] kmer = new B3Nucleotide[k];
		for(int i=0; i<k; i++)
			kmer[i] = symbols[its[i]];
		return kmer;
	}

	@Override
	public void kmer(B3Nucleotide[] ns) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int multiplicity() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public int[] positions() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public int[] sortedPositions() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public boolean isMinimalHapax() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public boolean isGlobalMaximalRepeat() {
		throw new UnsupportedOperationException("Unsupported function.");
	}

	@Override
	public boolean next() {
//		System.out.println("S "+it+" "+its[it]);
		
		B3Nucleotide p = B3Nucleotide.N;
		if(its[it] > -1){
			choosen[its[it]] = false;
			p = symbols[its[it]];
//			System.out.println(" p "+p);
		}
		
		its[it]++;
		while(its[it]<symbols.length && (choosen[its[it]] || symbols[its[it]]==p))
			its[it]++;
		
		if(its[it] == symbols.length){
			if(it == 0)
				return false;
			else{
				its[it] = -1;
				it--;
				return next();
			}
		}
		else{
			if(it == its.length-1){
				return true;
			}
			else{
				choosen[its[it]] = true;
				it++;
				return next();
			}
		}
		
	}

	@Override
	public boolean prev() {
		throw new UnsupportedOperationException("Not yet supported.");
	}

	@Override
	public boolean good() {
		throw new UnsupportedOperationException("Not yet supported.");
	}

	@Override
	public boolean hasNext() {
		throw new UnsupportedOperationException("Not yet supported.");
	}

	@Override
	public boolean hasPrev() {
		throw new UnsupportedOperationException("Not yet supported.");
	}

	@Override
	public int compare(IELSAIterator it) {
		throw new UnsupportedOperationException("Not yet supported.");
	}

	@Override
	public IELSAIterator clone(){
		throw new UnsupportedOperationException("Not yet supported.");
	}
}

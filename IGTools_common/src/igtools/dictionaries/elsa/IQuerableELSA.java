package igtools.dictionaries.elsa;

import igtools.common.sequence.B3Sequence;

public interface IQuerableELSA{

	public boolean contains(B3Sequence kmer);
	public IELSAIterator find(B3Sequence kmer);
	public IELSAIterator findOrNear(B3Sequence kmer, boolean perferForward);
	
	public int firstSAIndex(int pl, int pr);
	public int lastSAIndex(int pl, int pr);
	
	public int toUnique(int pl);
	public IELSAIterator toUniqueIterator(int pl);
	
	public int toUnique(int pl, int pr);
}

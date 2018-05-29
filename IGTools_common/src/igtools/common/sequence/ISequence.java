package igtools.common.sequence;

import igtools.common.kmer.IArrayKmer;

public interface ISequence <T extends IArrayKmer> {

	public int length();
	public T get(int pos, int length);
}

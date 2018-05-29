package igtools.common.streams;

import igtools.common.exceptions.StreamException;
import igtools.common.nucleotide.INucleotide;

public interface INucleotideStream<T extends INucleotide> {
	
	public static int NULL_CODE = 0xffffffff;
	

	public T next() 				throws StreamException;
	public int nextCode()				throws StreamException;
	
	public int position();
	public void seek(int position) 	throws StreamException;
	public int length();
	
	public void close()				throws StreamException;
	
}

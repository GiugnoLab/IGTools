package igtools.dictionaries.elsa;

import igtools.common.sequence.B3Sequence;

public interface IELSA {
	
	public int length();
	
	public B3Sequence b3seq();
	public int[] sa();
	public int[] lcp();
	
	
	public IELSAIterator begin(int k);
//	public IELSAIterator end(int k);
	
	public IELSAIterator rbegin(int k);
//	public IELSAIterator rend(int k);
	
	public int nof(int k);
	public int nof_mults(int k);
	
	public void load(String file) throws Exception;
	public void save(String file) throws Exception;
	
	public void print(int k);
	
}

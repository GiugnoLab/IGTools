package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;
import igtools.dictionaries.elsa.DiskNELSATreeNavigator;
import igtools.dictionaries.elsa.IELSAIterator;

public class DiskNELSATreeUnionNavigator {

	
	private DiskNELSATreeNavigator[] navigators = null;
	
	public DiskNELSATreeUnionNavigator(DiskNELSATreeNavigator[] navigators){
		this.navigators = navigators;
	}
	
	
	public void close() throws Exception{
		if(navigators != null){
			for(int i=0; i<navigators.length; i++){
				navigators[i].close();
			}
		}
	}
	
	
	public int nofNavigators(){
		if(navigators != null){
			return navigators.length;
		}
		return 0;
	}
	
	public long multiplicity(){
		long m = 0;
		for(int i=0; i<navigators.length; i++){
			if(navigators[i] != null){
				m += navigators[i].multiplicity();
			}
		}
		return m;
	}
	
	public void multiplicities(long[] m){
		for(int i=0; i<navigators.length; i++){
			if(navigators[i] != null){
				m[i] = navigators[i].multiplicity();
			}
			else{
				m[i] = 0;
			}
		}
	}
	public long[] multiplicities(){
		long[] m = new long[navigators.length];
		multiplicities(m);
		return m;
	}
	
	
	public DiskNELSATreeUnionNavigator getChild(int code) throws Exception{
		DiskNELSATreeNavigator[] c_navigators = new DiskNELSATreeNavigator[navigators.length];
		int count_null = 0;
		for(int i=0; i<navigators.length; i++){
			if(navigators[i] != null){
				c_navigators[i] = navigators[i].getChild(code);
				if(c_navigators[i] == null){
					count_null++;
				}
			}
			else{
				count_null++;
			}
		}
		if(count_null == navigators.length){
			return null;
		}
		else{
			return new DiskNELSATreeUnionNavigator(c_navigators);
		}
	}
	
	public static DiskNELSATreeUnionNavigator begin(B3Sequence[] b3seqs, String[] nelsas) throws Exception{
		DiskNELSATreeNavigator[] navigators = new DiskNELSATreeNavigator[b3seqs.length];
		for(int i=0; i<b3seqs.length; i++){
			navigators[i] = DiskNELSATreeNavigator.begin(b3seqs[i], nelsas[i]);
		}
		return new DiskNELSATreeUnionNavigator(navigators);
	}
}

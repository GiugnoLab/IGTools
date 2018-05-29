package igtools.analyses.recurrences.carpena;

import igtools.common.ds.StaticB2Tree;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;




/**
 * Modified Carpena2009 algorithm with parent comparison selection
 * 
 * @author vbonnici
 *
 */
public class Carpena2009ByParentC extends Carpena2009{

	public Carpena2009ByParentC(NELSA nelsa) {
		super(nelsa);
	}
	
	
	@Override
	public int select(int initial_k, double C0, boolean usePercentile, double C_delta, StaticB2Tree words){
		int count = 0;
		IELSAIterator  bit = nelsa.begin(initial_k);
		while(bit.next()){
			double C = C(bit);
			if((!usePercentile) || (C >= C0 - (C0 * C_delta))){
				count += recursive_select(initial_k, bit, C, C_delta, words);
			}
		}
		return count;
	}
	
	@Override
	protected int recursive_select(int initial_k, IELSAIterator it, double C0, double C_delta, StaticB2Tree words){
		int count = 0;
		double[] cs = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
		ExtensionNELSAIterator[] cs_its = {null, null, null, null};
		int i;
		double max = Double.NEGATIVE_INFINITY;
		
		
//		if(it!=null)
//			System.out.println(B3Nucleotide.toString(it.kmer()));
//		else
//			System.out.println(it);
		
		ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
		i=0;
		while(eit.next()){
			cs_its[i] = (ExtensionNELSAIterator)eit.clone();
			cs[i] = C(cs_its[i]);
//			System.out.println("cs[i] "+cs[i]+" > max " +max);
			if(cs[i] > max){
//				System.out.println("!");
				max = cs[i];
			}
			i++;
		}
		
//		System.out.println("max "+max);
//		System.out.print("C: ");
//		for(i=0; i<cs.length; i++){
//			if(cs[i] != Double.NEGATIVE_INFINITY){
//				System.out.print(i+"."+B3Nucleotide.toString(cs_its[i].kmer())+" ");
//			}
//		}
//		System.out.println();
//		System.out.print("C: ");
//		for(i=0; i<cs.length; i++){
//			if(cs[i] != Double.NEGATIVE_INFINITY){
//				System.out.print(cs[i]+" ");
//			}
//			if(cs[i] == Double.NaN)
//				System.exit(1);
//		}
//		System.out.println();
//		System.out.println();
		
//		System.out.println(max +" >= " + (C0 - (C0*C_delta)));
		if(i>0 &&      max >= (C0 - (C0*C_delta)) ){
			for(i=0; i<cs.length; i++){
				//System.out.println(cs[i]);
				if(cs[i] >= max){
					//System.out.println("select child "+i);
					if(cs_its[i].multiplicity()==1){
						words.add(nelsa.b3seq(), sa[cs_its[i].istart()], sa[cs_its[i].istart()] + it.k()+1);
						count++;
					}
					else{
						count += recursive_select(initial_k, cs_its[i], cs[i], C_delta, words);
					}
				}
			}
		}
		else if(it.k() != initial_k){
			words.add(nelsa.b3seq(), sa[it.istart()], sa[it.istart()] + it.k());
			count++;
		}
		return count;
		
	}
}

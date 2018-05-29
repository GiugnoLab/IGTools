package igtools.analyses.elongation;

import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class L2RElongation extends ElongationAlgorithm{

	public L2RElongation(NELSA nelsa) {
		super(nelsa);
	}

	@Override
	public void elongate(int startK, IteratorMeasure im, ElongationListener el) {
		IELSAIterator it = this.nelsa.begin(startK);
		Double it_value;
		while(it.next()){
			it_value = im.value(it);
			if(it_value != Double.NaN)
				elongate(it, it_value, im, el);
		}
	}
	
	private void elongate(IELSAIterator it, Double it_value, IteratorMeasure im, ElongationListener el) {
		
//		System.out.println("# "+ B3Nucleotide.toString(it.kmer())+" "+it_value);
		
		if(it_value != Double.NaN){
			ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
			Double eit_value;
			boolean notElongated = true;
			while(eit.next()){
				eit_value = im.value(eit);
				
//				System.out.println("> "+ B3Nucleotide.toString(eit.kmer())+" "+eit_value);
				
				if(eit_value > it_value){
					el.event(eit, ElongationEvent.EXTENDED, it_value, eit_value);
					notElongated = false;
					elongate(eit.clone(), eit_value, im, el);
				}
			}
			if(notElongated){
				el.event(it, ElongationEvent.MAXIMAL, it_value, it_value);
			}
		}
	}

}

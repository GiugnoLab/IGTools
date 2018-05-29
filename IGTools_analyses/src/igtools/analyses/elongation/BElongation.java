package igtools.analyses.elongation;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class BElongation extends ElongationAlgorithm{

	public BElongation(NELSA nelsa) {
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
	
	private static B3Nucleotide[] ee = {B3Nucleotide.A, B3Nucleotide.C, B3Nucleotide.G, B3Nucleotide.T};
	
	private void elongate(IELSAIterator it, Double it_value, IteratorMeasure im, ElongationListener el) {
		
		System.out.println("# "+ B3Nucleotide.toString(it.kmer())+" "+it_value);
		
		if(it_value != Double.NaN){
			
			B3Nucleotide[] ext = new B3Nucleotide[it.k()+2];
			it.kmer(ext);
			for(int i=ext.length-2; i>0; i--){
				ext[i] = ext[i-1];
			}
			
			boolean notElongated = true;
			IELSAIterator eit;
			Double eit_value;
			
			for(B3Nucleotide l : ee){
				ext[0] = l;
				for(B3Nucleotide r : ee){
					ext[ext.length-1] = r;
					eit = nelsa.find(ext);
					if(eit != null){
						eit_value = im.value(eit);
						
						System.out.println("> "+ B3Nucleotide.toString(ext)+" "+eit_value);
						if(eit_value > it_value){
							el.event(eit, ElongationEvent.EXTENDED, it_value, eit_value);
							notElongated = false;
							elongate(eit.clone(), eit_value, im, el);
						}
						
					}
//					else{
//						System.out.println("N: "+B3Nucleotide.toString(ext));
//					}
				}
			}
			
			if(notElongated){
				ext = new B3Nucleotide[it.k()+1];
				
				//L
				it.kmer(ext);
				for(int i=ext.length-1; i>0; i--){
					ext[i] = ext[i-1];
				}
				for(B3Nucleotide l : ee){
					ext[0] = l;
					eit = nelsa.find(ext);
					if(eit != null){
						eit_value = im.value(eit);
						
						System.out.println("> "+ B3Nucleotide.toString(ext)+" "+eit_value);
						if(eit_value > it_value){
							el.event(eit, ElongationEvent.EXTENDED, it_value, eit_value);
							notElongated = false;
							elongate(eit.clone(), eit_value, im, el);
						}
						
					}
//					else{
//						System.out.println("N: "+B3Nucleotide.toString(ext));
//					}
				}
				
				//R
				it.kmer(ext);
				for(B3Nucleotide r : ee){
					ext[ext.length-1] = r;
					eit = nelsa.find(ext);
					if(eit != null){
						eit_value = im.value(eit);
						
						System.out.println("> "+ B3Nucleotide.toString(ext)+" "+eit_value);
						if(eit_value > it_value){
							el.event(eit, ElongationEvent.EXTENDED, it_value, eit_value);
							notElongated = false;
							elongate(eit.clone(), eit_value, im, el);
						}
						
					}
					else{
						System.out.println("N: "+B3Nucleotide.toString(ext));
					}
				}
				
			}
			
			if(notElongated){
				el.event(it, ElongationEvent.MAXIMAL, it_value, it_value);
			}
		}
	}
}

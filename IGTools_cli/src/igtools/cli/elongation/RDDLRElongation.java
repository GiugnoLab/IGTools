package igtools.cli.elongation;

import igtools.analyses.elongation.BElongation;
import igtools.analyses.elongation.ElongationAlgorithm;
import igtools.analyses.elongation.L2RElongation;
import igtools.analyses.elongation.LnRElongation;
import igtools.analyses.elongation.R2LElongation;
import igtools.analyses.elongation.ElongationAlgorithm.ElongationEvent;
import igtools.analyses.recurrences.expcomp.KExpComp;
import igtools.common.distributions.distance.DistributionDistance;
import igtools.common.distributions.distance.KLDistance;
import igtools.common.distributions.distance.KSDistance;
import igtools.common.distributions.distance.PointDiffSUMDistance;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;


public class RDDLRElongation {
	public static void usage(){
		System.out.println("Usage: cmd b3seq nelsa <kl|ks|sum|ekl> <L2R|R2L|LnR|B> startK minMult ");
	}
	
	
	enum DistType {KL,KS,SUM, EKL};
	enum ElongType {L2R, R2L,LnR, B};
	
	public static void main(String[] args){
		String a_b3seq = "";
		String a_nelsa = "";
		DistType a_dtype = DistType.KL; 
		ElongType a_etype = ElongType.L2R;
		int k = 0;
		int a_minmult = 0;
		
		try{
			a_b3seq = args[0];
			a_nelsa = args[1];
			
			
			if(args[2].compareTo("kl") == 0){
				a_dtype = DistType.KL;
			}
			else if(args[2].compareTo("ks") == 0){
				a_dtype = DistType.KS;
			} 
			else if(args[2].compareTo("sum") == 0){
				a_dtype = DistType.SUM;
			}
			else if(args[2].compareTo("ekl") == 0){
				a_dtype = DistType.EKL;
			}
			else
				throw new Exception();
			
			
			if(args[3].compareTo("L2R")==0){
				a_etype = ElongType.L2R;
			}
			else if(args[3].compareTo("R2L")==0){
				a_etype = ElongType.R2L;
			}
			else if(args[3].compareTo("LnR")==0){
				a_etype = ElongType.LnR;
			}
			else if(args[3].compareTo("B")==0){
				a_etype = ElongType.B;
			}
			else
				throw new Exception();
			
			
			k = Integer.parseInt(args[4]);
			if(k<1)
				throw new Exception();
			
			a_minmult = Integer.parseInt(args[5]);
			if(k<0)
				throw new Exception();
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_b3seq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_nelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			DistributionDistance dd;
			switch(a_dtype){
			case KL:
				dd = new KLDistance.MaxKLDistance();
				break;
			case KS:
				dd = new KSDistance();
				break;
			case SUM:
				dd = new PointDiffSUMDistance();
				break;
			case EKL:
				dd = new KLDistance();
				break;
			default:
				dd = null;
				break;
			}
			
			final DistributionDistance fdd = dd;
			final int fstartK = k;
			final int minMult = a_minmult;
			
			
			ElongationAlgorithm.IteratorMeasure im = null;
			if(a_dtype == DistType.EKL){
				im = new ElongationAlgorithm.IteratorMeasure() {
					@Override
					public Double value(IELSAIterator it) {
						if(it.multiplicity() >= minMult)
							return KExpComp.exponentialDistance(it, fdd, 0.01);
						else
							return 0.0;
					}
				};
			}
			else{
				im = new ElongationAlgorithm.IteratorMeasure() {
					@Override
					public Double value(IELSAIterator it) {
						if(it.multiplicity() >= minMult)
							return KExpComp.distanceToExponential(it, fdd, 0.01);
						else
							return 0.0;
					}
				};
			}
			
			
			final ElongationAlgorithm.ElongationListener el = 
			new ElongationAlgorithm.ElongationListener() {
				@Override
				public void event(IELSAIterator it, ElongationEvent ev, Double parentValue, Double itValue) {
					if(ev == ElongationEvent.EXTENDED){
//						System.out.println("E: "+B3Nucleotide.toString(it.kmer())+" "+parentValue+" "+itValue);
					}
					else{
						if(it.k() != fstartK){
							System.out.println("S: "+B3Nucleotide.toString(it.kmer())+" "+itValue);
						}
					}
				}
					
			};
			
			
			ElongationAlgorithm ealgo = null;
			switch(a_etype){
			case B:
				ealgo = new BElongation(nelsa);
				break;
			case L2R:
				ealgo = new L2RElongation(nelsa);
				break;
			case LnR:
				ealgo = new LnRElongation(nelsa);
				break;
			case R2L:
				ealgo = new R2LElongation(nelsa);
				break;
			default:
				break;
			}
			
			ealgo.elongate(k, im, el);
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}

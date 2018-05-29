package igtools.cli.dictionaries;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class EnumerateKmers {

	
	public static void usage(){
		System.out.println("Usage: cmd b3seq nelsa k");
	}
	
	
	public static void main(String[] args){
		String a_b3seq = "";
		String a_nelsa = "";
		int k = 0;
		
		try{
			a_b3seq = args[0];
			a_nelsa = args[1];
			
			k = Integer.parseInt(args[2]);
			if(k<1)
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
			
			
//			IELSAIterator itf;
//			IELSAIterator it = nelsa.begin(k);
//			while(it.next()){
//				itf = nelsa.find(new B3LLSequence(it.kmer()));
//				System.out.println(it.istart()+"\t"+it.iend()+"\t"+B3Nucleotide.toString(it.kmer()));
//				System.out.println(itf.istart()+"\t"+itf.iend()+"\t"+B3Nucleotide.toString(itf.kmer()));
//			}
			B3Nucleotide[] kmer = new B3Nucleotide[k];
			IELSAIterator it = nelsa.begin(k);
			while(it.next()){
				it.kmer(kmer);
				System.out.println("K: "+B3Nucleotide.toString(kmer)+" "+it.multiplicity());
			}
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}

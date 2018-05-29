package igtools.cli.entropy;

import java.io.File;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import igtools.analyses.Inform;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;
import igtools.dictionaries.intersection.BigNumericUnion;
import igtools.dictionaries.intersection.NumericUnion;

public class DiskKEntropies {

	
	
	public static void usage(){
		System.out.println("Usage: csv[iseq.3bit iseq.nelsa] [ks]");
	}
	
	public static void main(String[] args){

		try{
			String a_config = args[0];
			
			Vector<String> iseqs = new Vector<String>();
			Vector<String> inelsas = new Vector<String>();
			int nof_chrs = 0;
			String scanner_line;
			String[] scanner_cols;
			Scanner scanner = new Scanner(new File(a_config));
			while(scanner.hasNextLine()){
				scanner_line = scanner.nextLine().trim();
				if(scanner_line.length() > 0){
					nof_chrs++;
					scanner_cols = scanner_line.split(" ");
					iseqs.add(scanner_cols[0]);
					inelsas.add(scanner_cols[1]);
				}
			}
			scanner.close();
			
			Set<Integer> ks = new TreeSet<Integer>();
			for(int i=1; i<args.length; i++){
				if(args[i].contains("-")){
					int k1 = Integer.parseInt(args[i].split("-")[0]);
					int k2 = Integer.parseInt(args[i].split("-")[1]);
					for(int j=k1; j<=k2; j++){
						ks.add(j);
					}
				}
				else{
					ks.add(Integer.parseInt(args[i]));
				}
			}
			
			B3LLSequence[] seqs = new B3LLSequence[nof_chrs];
			
			for(int i=0; i<nof_chrs; i++){
				System.out.println(iseqs.get(i));			
				seqs[i] = B3LLSequence.load(iseqs.get(i));
			}
			
			double length = 0.0;
			double bads = 0.0;
			double a_length = 0.0;
			
			for(int i=0; i<nof_chrs; i++){
				System.out.println(i+" "+iseqs.get(i)+" "+seqs[i].length()+" "+seqs[i].countBads());
				length += seqs[i].length();
				bads += seqs[i].countBads();
			}
			a_length = length - bads;
			
			double log_2 = Math.log(length) / Math.log(2);
			double m = 2.0 * Math.log(length) / Math.log(4);
			double a_log_2 = Math.log(a_length) / Math.log(2);
			double a_m = 2.0 * Math.log(a_length) / Math.log(4);
			
			System.out.println("m "+m);
			System.out.println("a_m "+a_m);
			System.out.println("log_2 "+log_2);
			
			
//			System.out.println("auto added: "+((int)Math.floor(m))+" "+((int)Math.floor(m) - 1)+" "+((int)Math.ceil(m))+" "+((int)Math.ceil(m) + 1));
//			ks.add((int)Math.ceil(m));
//			ks.add((int)Math.floor(m));
//			ks.add((int)Math.ceil(m) + 1);
//			ks.add((int)Math.floor(m) - 1);
			
			for(int k : ks){
				double[] hs = entropy(nof_chrs, inelsas, seqs, k, length);
				System.out.println("E_t: "+k+" "+hs[0]);
				System.out.println("E_b: "+k+" "+hs[1]);
			}
//			for(int k=min_k; k<= (int)Math.ceil(m) + 1; k++){
//				double h = entropy(nof_chrs, inelsas, seqs, k);
//				System.out.println("E: "+k+" "+h+" "+(h-log_2));
//				//double h_m_ceil = entropy(nof_chrs, inelsas, seqs, (int)Math.ceil(m));
//				//System.out.println("bm_ceil "+ (h_m_ceil - log_2));
//				//double h_m_floor = entropy(nof_chrs, inelsas, seqs, (int)Math.floor(m));
//				//System.out.println("bm_floor "+ (h_m_floor - log_2));
//			}
			
		}catch(Exception e){
			e.printStackTrace(System.err);
			System.err.println(e);
		}
		
	}
	
	public static double[] entropy(int nof_chrs, Vector<String> inelsas, B3LLSequence[] seqs, int k, double length) throws Exception{
		
		IELSAIterator[] its = new IELSAIterator[nof_chrs];
		for(int i=0; i<nof_chrs; i++){
			its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		}
		
		final long[] nof = {0L,0L,0L};
		final double[] count = {0.0, length -k + 1.0};
		
		final double[] e = {0.0,0.0};
		final double[] p = {0.0};
		
		
		
		NumericUnion.UnionListerner listener = new NumericUnion.UnionListerner() {
			@Override
			public void intersection(B3Nucleotide[] kmer) {
			}
			@Override
			public void union(B3Nucleotide[] kmer, int mult) {
				count[0] += (double)mult;
				nof[0]++;
				if(mult == 1){
					nof[1]++;
				}
				else if(mult > 1){
					nof[2]++;
				}
			}
		};
		
		NumericUnion union = null;
		if(k < 32)
			union = new NumericUnion(seqs, its, k, listener);
		else 
			union = new BigNumericUnion(seqs, its, k, listener);
		union.union();
		
		System.out.println("count "+count[0]);
		System.out.println("nof s "+nof[0]+" "+nof[1]+" "+nof[2]);
		
		for(int i=0; i<its.length; i++){
			if(its[i] != null)
			((OnDiskNELSAIteratorV2)its[i]).close();
			its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		}
		
		listener = new NumericUnion.UnionListerner() {
			@Override
			public void intersection(B3Nucleotide[] kmer) {
			}
			@Override
			public void union(B3Nucleotide[] kmer, int mult) {
				p[0] = ((double)mult) / count[0];
				e[0] += p[0] * Math.log(p[0]);
				
				p[0] = ((double)mult) / count[1];
				e[1] += p[0] * Math.log(p[0]);
			}
		};
		union = null;
		if(k < 32)
			union = new NumericUnion(seqs, its, k, listener);
		else 
			union = new BigNumericUnion(seqs, its, k, listener);
		union.union();
		
		
		for(int i=0; i<its.length; i++){
			if(its[i] != null)
			((OnDiskNELSAIteratorV2)its[i]).close();
		}
		
		for(int i=0; i<e.length; i++){
			e[i] *= -1;
			e[i] /= Math.log(2);
		}
		
//		double entr = e[0];
//		System.out.println("entr "+entr);
//		entr *= -1;
//		entr /= Math.log(2);
//		return entr;
		return e;
	}
}

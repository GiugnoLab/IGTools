package igtools.cli.wordset;

import igtools.analyses.sequence.coverage.DictionarySequenceCoverage;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * Calculate the total coverage of a set of sequences [seqs] (in FASTA format) in a given string [b3seq, nelsa] G.
 *   
 * @author vbonnici
 *
 */
public class Coverage {

	public static void usage(){
		System.out.println("IGTools: Calculate the total coverage of a set of sequences [seqs] (in FASTA format) in a given string [b3seq, nelsa]");
		System.out.println("Usage: cmd seqs b3seq nelsa");
		System.out.println("\t INPUT: seqs must be a list of sequences in FASTA format, one per line.");
	}
	
	public static void main(String[] args){
		String a_seqs = "";
		String a_b3seq = "";
		String a_nelsa = "";
		
		try{
			a_seqs = args[0];
			a_b3seq = args[1];
			a_nelsa = args[2];
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
			
			
			//boolean[] coverage = new boolean[b3seq.length()];
			//Arrays.fill(coverage, false);
			double[] coverage = new double[b3seq.length()];
			Arrays.fill(coverage,0);
			
			
			IELSAIterator it;
			int[] sa = nelsa.sa();
			int k;
			
			System.out.println("reading...");

			Scanner scanner = new Scanner(new File(a_seqs));
			String line;
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					it = nelsa.find(new B3LLSequence(line));
					
					if(it != null){
						k = it.k();
						for(int i= it.istart(); i<it.iend(); i++){
							for(int j=sa[i]; j<sa[i]+k; j++){
								//coverage[j] = true;
								coverage[j]++;
							}
						}
					}
				}
			}
			scanner.close();
			
			System.out.println("calculating...");
			double tot = (double)(b3seq.length() - b3seq.countBads());
			double cov = 0.0;
			for(int i=0; i<coverage.length; i++){
				//if(coverage[i] == true){
				if(coverage[i] > 0){
					cov++;
				}
			}
			
			double avg_cov = igtools.common.util.AvgSD.avg(coverage, 0.0);
			double sd_cov = igtools.common.util.AvgSD.sd(coverage, avg_cov, 0.0);
			
			System.out.println("coverage "+(cov));
			System.out.println("coverage ratio "+(cov/tot));
			System.out.println("avg coverage "+avg_cov);
			System.out.println("sd coverage "+sd_cov);
			
			Map<Double, Double> covDistr = new TreeMap<Double,Double>();
			DictionarySequenceCoverage.coverageDistribution(coverage, covDistr);
			
			for(Map.Entry<Double, Double> en : covDistr.entrySet()){
				System.out.println("CD "+en.getKey()+" "+en.getValue());
			}
			
			System.out.println(".");
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}

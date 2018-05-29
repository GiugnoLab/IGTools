package igtools.cli.wordset;

import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

public class CoverageByLength {

	public static void usage(){
		System.out.println("Usage: cmd seqs b3seq nelsa");
		System.out.println("\t INPUT: seqs must be a list of sequences in FASTA format, one per line.");
		System.out.println("\t OUTPUT: a file with <LCD: length count seqCov covpAvg covpSD> lines.");
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
			
			
			Map<Integer,Double> m_count = new TreeMap<Integer,Double>();
			
			
			
			IELSAIterator it;			
			String line;
			
			System.out.println("reading...");
			Scanner scanner = new Scanner(new File(a_seqs));
			while(scanner.hasNextLine()){
				try{
					line = scanner.nextLine().trim();
					if(line.length() > 0){
						addTo(line.length(), 1.0, m_count);
						//it = nelsa.find(new B3LLSequence(line));
					}
				}catch(Exception e){
					e.printStackTrace();
					System.out.println(e);
				}
			}
			scanner.close();
			
			
			double[] coverage = new double[b3seq.length()];
			Arrays.fill(coverage,0);
			
			Double c;
			int k;
			int[] sa = nelsa.sa();
			
			double scov, covpAvg, covpSD;
			double tot = (double)(b3seq.length() - b3seq.countBads());
			
			for(Integer length : m_count.keySet()){
				System.out.println("L "+length);
				c = m_count.get(length);
				Arrays.fill(coverage,0);
				
				scanner = new Scanner(new File(a_seqs));
				while(scanner.hasNextLine()){
					try{
						line = scanner.nextLine().trim();
						if(line.length() > 0 && line.length()==length){
							
							it = nelsa.find(new B3LLSequence(line));
							
							if(it != null){
								k = it.k();
								
								if(it != null){
									for(int i= it.istart(); i<it.iend(); i++){
										for(int j=sa[i]; j<sa[i]+k; j++){
											//coverage[j] = true;
											coverage[j]++;
										}
									}
								}
							}
						}
					}catch(Exception e){
						e.printStackTrace();
						System.out.println(e);
					}
				}
				scanner.close();
				
				scov = 0.0;
				for(int i=0; i<coverage.length; i++){
					if(coverage[i] > 0){
						scov++;
					}
				}
				scov /= tot;
				
				covpAvg = igtools.common.util.AvgSD.avg(coverage, 0.0);
				covpSD = igtools.common.util.AvgSD.sd(coverage, covpAvg, 0.0);
				
				System.out.println("LCD: "+length+" "+c+" "+scov+" "+covpAvg+" "+covpSD);
				
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
	
	private static void addTo(Integer key, Double value, Map<Integer,Double> d){
		Double v = d.get(key);
		if(v == null){
			v = 0.0;
		}
		v += value;
		d.put(key, v);
	}
}

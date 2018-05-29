package igtools.cli.entropy;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.CompleteIterator;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;
import igtools.dictionaries.intersection.DiskNELSATreeUnionNavigator;
import igtools.dictionaries.intersection.NumericUnion;

import java.io.File;
import java.util.Scanner;
import java.util.Vector;

public class DiskJKL {

	
	public static void usage(){
		System.out.println("Usage: cmd csv[seq nelsa] csv[seq nelsa]");
	}
	
	public static void main(String[] args){
		try{
			long start_t;
			
			Vector<String> iseqs1 = new Vector<String>();
			Vector<String> inelsas1 = new Vector<String>();
			Vector<String> iseqs2 = new Vector<String>();
			Vector<String> inelsas2 = new Vector<String>();
			getConf(args[0], iseqs1, inelsas1);
			getConf(args[1], iseqs2, inelsas2);
			
			System.out.println("Loading seqs 1...");
start_t=System.currentTimeMillis();
			B3LLSequence[] seqs1 = getSeqs(iseqs1);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			System.out.println("Loading seqs 2...");
start_t=System.currentTimeMillis();
			B3LLSequence[] seqs2 = getSeqs(iseqs2);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			
			String[] nelsas1 = new String[inelsas1.size()]; inelsas1.toArray(nelsas1);
			String[] nelsas2 = new String[inelsas2.size()]; inelsas2.toArray(nelsas2);
			
			
			int m1 =  (int)Math.ceil(Math.log(totalLength(seqs1)) / Math.log(2));
			int m2 =  (int)Math.ceil(Math.log(totalLength(seqs2)) / Math.log(2));
			int m = 0;
			if(m1 <= m2) m = m1;
			else m = m2;
			
			System.out.println("ceil(log2("+totalLength(seqs1)+")) = "+m1+", ceil(log2("+totalLength(seqs2)+")) = "+m2+"");
			System.out.println("m "+ m);
			
			final double[] nofs = {0.0};
			final double[] sums = {0.0,0.0};
			System.out.println("Counts...");
			//ms_counts(nav1, nav2, nofs, sums, 0, m);
			final double[] sum = {0.0};
			System.out.println("seqs 1...");
start_t=System.currentTimeMillis();
			its_counts(seqs1, nelsas1, m, nofs, sum);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			System.out.println(nofs[0]+" "+sum[0]);
			sums[0] = sum[0];
			sum[0] = 0.0;
			System.out.println("seqs 2...");
start_t=System.currentTimeMillis();
			its_counts(seqs2, nelsas2, m, nofs, sum);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			System.out.println(nofs[0]+" "+sum[0]);
			sums[1] = sum[0];
			
			
			
			
			System.out.println("Init navigator 1...");
start_t=System.currentTimeMillis();
			DiskNELSATreeUnionNavigator nav1 = DiskNELSATreeUnionNavigator.begin(seqs1, nelsas1);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			System.out.println("Init navigator 2...");
start_t=System.currentTimeMillis();
			DiskNELSATreeUnionNavigator nav2 = DiskNELSATreeUnionNavigator.begin(seqs2, nelsas2);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);
			System.out.println("JKL...");			
			double[] jkl = {0.0};
start_t=System.currentTimeMillis();
			jkl(nav1, nav2, sums, 0, m, jkl);
System.out.println("time "+(System.currentTimeMillis()-start_t)/1000);			
			
			System.out.println("JKL "+jkl[0]);
			
			nav1.close();
			nav2.close();
			
		}catch(Exception e){
			e.printStackTrace(System.err);
			System.err.println(e);
		}
	}
	
	
	public static void its_counts(B3LLSequence[] seqs, String[] nelsas, int k, final double[] count, final double[] sum) throws Exception{
		OnDiskNELSAIteratorV2[] its = diskIterators(seqs, nelsas, k);
		NumericUnion.UnionListerner listener = new NumericUnion.UnionListerner() {
			@Override
			public void intersection(B3Nucleotide[] kmer) {
			}
			@Override
			public void union(B3Nucleotide[] kmer, int mult) {
//				System.out.println(B3Nucleotide.toString(kmer));
				sum[0] += (double)mult;
				count[0]++;
			}
		};
		NumericUnion union = new NumericUnion(seqs, its, k, listener);
		union.union();
	}
	
	
	public static void ms_counts(DiskNELSATreeUnionNavigator nav1, DiskNELSATreeUnionNavigator nav2, double[] nofs, double[] sums, int k, int m)
			throws Exception
	{
		if(k == m){
			nofs[0]++;
			sums[0] += nav1.multiplicity();
			sums[1] += nav2.multiplicity();
			return;
		}
		else{
			DiskNELSATreeUnionNavigator[] c1 = new DiskNELSATreeUnionNavigator[4];
			DiskNELSATreeUnionNavigator[] c2 = new DiskNELSATreeUnionNavigator[4];
			
			for(int c=0; c<4; c++){
				c1[c] = nav1.getChild(c);
				c2[c] = nav2.getChild(c);
				if((c1[c]==null && c2[c]!=null) || (c1[c]!=null && c2[c]==null)){
					nofs[0]++;
					sums[0] += nav1.multiplicity();
					sums[1] += nav2.multiplicity();
					return;
				}
			}
			
			for(int c=0; c<4; c++){
				if(c1[c] != null){
					ms_counts(c1[c], c2[c], nofs, sums, k+1, m);
				}
			}
		}
	}
	
	
	private static final double log2 = Math.log(2.0);
	public static void jkl(DiskNELSATreeUnionNavigator nav1, DiskNELSATreeUnionNavigator nav2, final double[] sums, int k, int m, double[] jj)
			throws Exception
	{
		if(k == m){
			double p0 = ((double)nav1.multiplicity()) / sums[0];
			double p1 = ((double)nav2.multiplicity()) / sums[1];
			jj[0] += (p0 - p1) * Math.log(p0 / p1) * log2; 
			return;
		}
		else{
			DiskNELSATreeUnionNavigator[] c1 = new DiskNELSATreeUnionNavigator[4];
			DiskNELSATreeUnionNavigator[] c2 = new DiskNELSATreeUnionNavigator[4];
			
			for(int c=0; c<4; c++){
				c1[c] = nav1.getChild(c);
				c2[c] = nav2.getChild(c);
				if((c1[c]==null && c2[c]!=null) || (c1[c]!=null && c2[c]==null)){
					double p0 = ((double)nav1.multiplicity()) / sums[0];
					double p1 = ((double)nav2.multiplicity()) / sums[1];
					jj[0] += (p0 - p1) * Math.log(p0 / p1) * log2;
					return;
				}
			}
			
			for(int c=0; c<4; c++){
				if(c1[c] != null){
					jkl(c1[c], c2[c], sums, k+1, m, jj);
				}
			}
		}
	}
	
	
	
	
	private static long totalLength(B3LLSequence[] seqs){
		long l = 0l;
		for(int i=0; i<seqs.length; i++){
			l += seqs[i].length();
		}
		return l;
	}
	private static void getConf(String file, Vector<String> seqs, Vector<String> nelsas) throws Exception{
		int nof_chrs = 0;
		String scanner_line;
		String[] scanner_cols;
		Scanner scanner = new Scanner(new File(file));
		while(scanner.hasNextLine()){
			scanner_line = scanner.nextLine().trim();
			if(scanner_line.length() > 0){
				nof_chrs++;
				scanner_cols = scanner_line.split(" ");
				seqs.add(scanner_cols[0]);
				nelsas.add(scanner_cols[1]);
			}
		}
		scanner.close();
	}
	private static B3LLSequence[] getSeqs(Vector<String> iseqs) throws Exception{
		B3LLSequence[] seqs = new B3LLSequence[iseqs.size()];
		for(int i=0; i< iseqs.size(); i++){
			System.out.println(iseqs.get(i));
			seqs[i] = B3LLSequence.load(iseqs.get(i));
		}
		return seqs;
	}
	
	private static OnDiskNELSAIteratorV2[] diskIterators(B3LLSequence[] seqs, String[] nelsas, int k) throws Exception{
		OnDiskNELSAIteratorV2[] its = new OnDiskNELSAIteratorV2[nelsas.length];
		for(int i=0; i<nelsas.length; i++){
			its[i] = new OnDiskNELSAIteratorV2(seqs[i], nelsas[i], k);
		}
		return its;
	}
	
}

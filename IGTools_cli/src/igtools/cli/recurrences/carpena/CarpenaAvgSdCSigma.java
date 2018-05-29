package igtools.cli.recurrences.carpena;

import igtools.analyses.recurrences.carpena.Carpena2009Factory;
import igtools.common.ds.StaticB2Tree;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;

import java.util.Map;
import java.util.TreeMap;

public class CarpenaAvgSdCSigma {
	private static void usage(){
		System.out.println("Usage: cmd sequence.3bit nelsa.nelsa [original|parent] k... ");
	}
	
	
	
	private static double avg(double a[]){
		double avg = 0.0;
		for(int i=0; i<a.length; i++){
			avg += a[i];
		}
		return avg / ((double)a.length);
	}
	private static double sd(double[] a, double avg){
		double sd = 0.0;
		for(int i=0; i<a.length; i++){
			sd += (a[i] - avg) * (a[i] - avg);
		}
		return Math.sqrt(sd / ((double) a.length));
	}
	
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_inelsa = "";
		String a_type = "";
		int[] ks = null;
		
		boolean original = false;		
		
		try{
			a_iseq = args[0];
			a_inelsa = args[1];
			
			a_type = args[2];
			if(a_type.compareTo("original")==0)
				original = true;
			else if(a_type.compareTo("parent")==0)
				original = false;
			else
				throw new Exception();
			
			if(args.length <= 3)
				throw new Exception();
			
			System.out.println(args.length +" "+ (args.length - 3));
			
			ks = new int[args.length - 3];
			for(int i=3; i<args.length; i++){
				ks[i-3] = Integer.parseInt(args[i]);
			}
		}catch(Exception e){
			e.printStackTrace();
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			StaticB2Tree words = new StaticB2Tree();
			//igtools.analyses.words.Carpena2009 c2009 = new igtools.analyses.words.Carpena2009(nelsa);
			igtools.analyses.recurrences.carpena.Carpena2009 c2009 = Carpena2009Factory.factory(nelsa, !original, false);
			
			System.out.println("k sigma_avg sigma_sd sigma_nor_avg sigma_nor_sd C_avg C_sd");
			
			for(int i=0; i<ks.length; i++){
				int k = ks[i];
				System.out.print(k+" ");
				
				double avg;
				double[] cs = c2009.get_sigmas(k);
				avg = avg(cs);
				System.out.print(avg + " " +sd(cs, avg) +" ");
				
				cs = c2009.get_sigma_nors(k);
				avg = avg(cs);
				System.out.print(avg + " " +sd(cs, avg) +" ");
				
				cs = c2009.get_Cs(k);
				avg = avg(cs);
				System.out.print(avg + " " +sd(cs, avg) +" ");
				
				System.out.println();
			}
			
		}catch(Exception ex){
			ex.printStackTrace();
		}
	}
}
package igtools.cli;

import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class GenomeKStats {

	
	public static void usage(){
		System.out.println("Usage: cmd from_k to_k iseq.3bit iseq.nelsa");
		System.out.println("for each value of k, it outputs the following informations separated by a space:");
		System.out.println("\tk |D_k| |H_k| |R_k| |T_k| |E_k|");
	}
	
	public static void main(String[] args){
		
		int from_k = 0;
		int to_k = 0;
		String a_iseq = "";
		String a_inelsa = "";
		
		try{
			from_k = Integer.parseInt(args[0]);
			if(from_k < 0) throw new Exception();
			to_k = Integer.parseInt(args[1]);
			if(to_k < 0 || to_k < from_k) throw new Exception();
			a_iseq = args[2];
			a_inelsa = args[3];
		}catch(Exception  e){
			usage();
			System.exit(0);
		}
		
		
		try{

			System.out.println("seq");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("nelsa");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			nelsa.setSequence(b3seq);
			

//			int k=1;
//			int max_k = 20;
//			while(true){
//				if(k>max_k)
//					break;
//				IELSAIterator it = nelsa.begin(k);
//				int tot=0,hapaxes=0, repeats=0;
//				if(it.next()){
//					do{
//						tot++;
//						if(it.multiplicity() == 1){
//							hapaxes++;
//						}
//						else{
//							repeats++;
//						}
//					}while(it.next());
//				}
//				else break;
//				
//				System.out.println(k+" "+Math.pow(4, k)+" "+( ((double)tot) /Math.pow(4, k))+ " "+tot+" "+hapaxes+" "+repeats+" "+(((double)hapaxes)/((double)tot)));
//				k++;
//			}
			
			System.out.println("k |D_k| |H_k| |R_k| |T_k| |E_k|");
			
			for(int k=from_k; k<=to_k; k++){
				int dk=0, hk=0, rk=0, tk=0;
				
				IELSAIterator it = nelsa.begin(k);
				while(it.next()){
					dk++;
					tk += it.multiplicity();
					if(it.multiplicity() == 1){
						hk++;
					}
					else{
						rk++;
					}
				}
				
				double dtk = (double)tk;
				double ek = 0.0;
				double p;
				it = nelsa.begin(k);
				while(it.next()){
					p = ((double)it.multiplicity()) / dtk;
					ek += p * Math.log(p);
				}
				ek /= Math.log(2);
				
				System.out.println(k+" "+dk+" "+hk+" "+rk+" "+tk+" "+ek);
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
}

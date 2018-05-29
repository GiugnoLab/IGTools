package igtools.cli.wordset;

import java.io.File;
import java.util.Scanner;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

/**
 * Cover by N the occurrences of the sequences in iseqs.fa
 * 
 * @author vbonnici
 *
 */
public class CoverByN {

	public static void usage(){
		System.out.println("Usage: cmd iseqs.fa iseq.3bit inelsa.nelsa oseq.3bit");
	}
	
	public static void main(String[] args){
		String a_fas = "";
		String a_iseq = "";
		String a_inelsa = "";
		String a_oseq = "";
		
		
		try{
			
			a_fas = args[0];
			a_iseq = args[1];
			a_inelsa = args[2];
			a_oseq = args[3];
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence iseq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			System.out.println("done");
			nelsa.setSequence(iseq);
			
			
			System.out.println("Cloning...");
			B3LLSequence oseq = iseq.clone();
			System.out.println("done");
			
			
			int ccount = 0;
			
			System.out.println("Covering...");
			IELSAIterator it;
			Scanner scanner = new Scanner(new File(a_fas));
			String line;
			int is,ie,il;
			int[] sa = nelsa.sa();
			while(scanner.hasNextLine()){
				try{
					line = scanner.nextLine().trim();
					if(line.length() > 0){
						il = line.length();
						it = nelsa.find(new B3LLSequence(line));
						if(it != null){
							is = it.istart();
							ie = it.iend();
							
							System.out.println(line+" "+il+" ["+is+","+ie+"]"+" "+(ie-is));
							
							for(int i=is; i<ie; i++){
								for(int j=0; j<il; j++){
									if(oseq.getB3(sa[i]+j)!=B3Nucleotide.N_CODE)
										ccount++;
									oseq.setB3(sa[i]+j, B3Nucleotide.N_CODE);
								}
							}
						}
						else{
							System.out.println("#");
						}
					}
				}catch(Exception e){
					e.printStackTrace();
					System.out.println(e);
				}
			}
			scanner.close();
			System.out.println("done");
			
			System.out.println("Saving...");
			oseq.save(a_oseq);
			System.out.println("done");
			
			System.out.println("s length "+iseq.length());
			System.out.println("s bads "+iseq.countBads());
			System.out.println("c count "+ccount);
			System.out.println("c/length "+( (double)ccount/ (double)iseq.length()));
			System.out.println("c/length-bads "+((double)ccount/ (double)(iseq.length()-iseq.countBads())));
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}

package igtools.cli.util;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.BufferedFASTANucleotideStream;
import igtools.common.streams.INucleotideStream;

public class Extract {


	private static void usage(){
		System.out.println("Usage: cmd fasta isequence.fasta start end osequence.fasta [--no-ns | --skip-ns | --ns-bridges]");
	}
	
	private enum NsPolicy{
		INLCUDE,
		NO_NS,
		SKIP_NS,
		NS_BRIDGES
	}
	
	public static void main(String[] args){
		
		
		String a_iseq = null;
		String a_oseq = null;
		NsPolicy a_npol = NsPolicy.INLCUDE;
		int a_start = -1;
		int a_end = -1;
		
		int line_length = 80;
		
		try{
			//args[0] = fasta
			
			a_iseq = args[1];
			a_start = Integer.parseInt(args[2]);
			a_end = Integer.parseInt(args[3]);
			if(a_end == -1)
				a_end = Integer.MAX_VALUE;
			
			a_oseq = args[4];
			
			if(args.length>5){
				if(args[5].compareTo("--no-ns")==0)
					a_npol = NsPolicy.NO_NS;
				else if(args[5].compareTo("--skip-ns")==0)
					a_npol = NsPolicy.SKIP_NS;
				else if(args[5].compareTo("--ns-bridges")==0)
					a_npol = NsPolicy.NS_BRIDGES;
				else{
					usage();
					System.exit(1);
				}
			}
			
		}catch(Exception e){
			usage();
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println("[ARGS][output file]["+a_oseq+"]");
		
		
		try{
			PrintWriter writer = new PrintWriter(new BufferedOutputStream(new FileOutputStream(a_oseq, false)));
			
			int i; 
			INucleotideStream<B3Nucleotide> fastaStream = new BufferedFASTANucleotideStream(a_iseq);
			int n = B3Nucleotide.NULL_CODE;
			
			System.out.println("Seeking...");
			i = 0;
			while(i!=a_start && fastaStream.nextCode() != INucleotideStream.NULL_CODE)
				i++;
			
			if(i!=a_start){
				System.out.println("No much symbols to start");
				System.exit(0);
			}
			
			
			System.out.println("Extracting...");
			
			
			if(a_npol == NsPolicy.NO_NS){
				System.out.println("NO Ns mode");
				int l = 0;
				while(i<=a_end && (n=fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
					if(n == B3Nucleotide.N_CODE){
						System.out.println("Found N at position "+i+" of the input sequence: stopping...");
						break;
					}
					writer.print(B3Nucleotide.charFor(n));
					i++;
					l++;
					
					if(l==line_length){
						writer.println("");
						l=0;
					}
				}
			}
			else if(a_npol == NsPolicy.NS_BRIDGES){
				boolean on_ns = false;
				int l = 0;
				while(i<=a_end && (n=fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
					if(on_ns){
						if(n != B3Nucleotide.N_CODE){
							on_ns = false;
							writer.print(B3Nucleotide.charFor(n));
							i++;
							l++;
							if(l==line_length){
								writer.println("");
								l=0;
							}
						}
					}
					else{
						if(n == B3Nucleotide.N_CODE)
							on_ns = true;
						writer.print(B3Nucleotide.charFor(n));
						i++;
						l++;
						if(l==line_length){
							writer.println("");
							l=0;
						}
					}
				}
			}
			else if(a_npol == NsPolicy.SKIP_NS){
				int l = 0;
				int skipped_ns = 0;
				while(i<=a_end && (n=fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
					if(n != B3Nucleotide.N_CODE){
						writer.print(B3Nucleotide.charFor(n));
						i++;
						l++;
						
						if(l==line_length){
							writer.println("");
							l=0;
						}
					}
					else{
						skipped_ns++;
					}
				}
				System.out.println("Skipped Ns: "+skipped_ns);
			}
			else{//INCLUDE
				int l = 0;
				while(i<=a_end && (n=fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
					writer.print(B3Nucleotide.charFor(n));
					i++;
					l++;
					
					if(l==line_length){
						writer.println("");
						l=0;
					}
				}
			}
			
			if(n == INucleotideStream.NULL_CODE){
				System.out.println("[Warning] No much symbols to end, out length: "+(i-a_start));
			}
			
			System.out.println("Total extract synbols: "+(i-a_start));
			
			writer.flush();
			writer.close();
			fastaStream.close();
			
		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
			System.exit(2);
		}
		System.exit(0);
	}
}

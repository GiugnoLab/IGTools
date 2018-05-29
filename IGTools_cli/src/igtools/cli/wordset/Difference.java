package igtools.cli.wordset;

import igtools.common.ds.StaticB2Tree;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.io.File;
import java.util.Scanner;

public class Difference {

	
	
	public static void usage(){
		System.out.println("Usage: cmd [exact|prefix|suffix|sub]  file1 file2");
		System.out.println("OUTPUT:  <D: word>");
	}
	
	public static void main(String[] args){
		
		String a_type = "";
		String a_file1 = "";
		String a_file2 = "";
		
		try{
			a_file1 = args[1];
			a_file2 = args[2];
			
			a_type = args[0];
			if(a_type.compareTo("exact")==0){
				exact(a_file1, a_file2);
			}
			else if(a_type.compareTo("prefix")==0){
				prefix(a_file1, a_file2);
			}
			else if(a_type.compareTo("suffix")==0){
				suffix(a_file1, a_file2);
			}
			else if(a_type.compareTo("sub")==0){
				sub(a_file1, a_file2);
			}
			else{
				throw new Exception();
			}
			
			
		}catch(Exception e){
			usage();
			System.exit(0);
		}
		
	}
	
	
	public static void exact(String file1, String file2){
		try{
			
			StaticB2Tree f2tree = new StaticB2Tree();
			
			String line;
			Scanner scanner = new Scanner(new File(file2));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					f2tree.add(new B3LLSequence(line), 0, line.length());
				}
			}
			scanner.close();
			
			scanner = new Scanner(new File(file1));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					if(!(f2tree.containsOriginal(new B3LLSequence(line)))){
						System.out.println("D: "+line);
					}
				}
			}
			scanner.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
	
	public static void prefix(String file1, String file2){
		try{
			
			StaticB2Tree f2tree = new StaticB2Tree();
			
			String line;
			Scanner scanner = new Scanner(new File(file2));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					f2tree.add(new B3LLSequence(line), 0, line.length());
				}
			}
			scanner.close();
			
			scanner = new Scanner(new File(file1));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					if(!(f2tree.contains(new B3LLSequence(line)))){
						System.out.println("D: "+line);
					}
				}
			}
			scanner.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
	
	
	public static void suffix(String file1, String file2){
		try{
			
			StaticB2Tree f2tree = new StaticB2Tree();
			
			String line;
			Scanner scanner = new Scanner(new File(file2));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					f2tree.addReverse(new B3LLSequence(line), 0, line.length());
				}
			}
			scanner.close();
			
			scanner = new Scanner(new File(file1));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					if(!(f2tree.contains(new B3LLSequence(line)))){
						System.out.println("D: "+line);
					}
				}
			}
			scanner.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
	
	public static void sub(String file1, String file2){
		try{
			
			String s2 = "";
			
			String line;
			Scanner scanner = new Scanner(new File(file2));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					s2 += line+"N";
				}
			}
			scanner.close();
			
			if(s2.length() > 0){
			
				B3LLSequence b2 = new B3LLSequence(s2);
				s2 = null;
				NELSA n2 = new NELSA(b2);
				
				IELSAIterator it;
				scanner = new Scanner(new File(file1));
				while(scanner.hasNextLine()){
					line = scanner.nextLine().trim();
					if(line.length() > 0){
						it = n2.find(new B3LLSequence(line));
						if(it == null){
							System.out.println("D: "+line);
						}
					}
				}
				scanner.close();
			}
			else{
				scanner = new Scanner(new File(file1));
				while(scanner.hasNextLine()){
					line = scanner.nextLine().trim();
					if(line.length() > 0){
						System.out.println("D: "+line);
					}
				}
				scanner.close();
			}
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}

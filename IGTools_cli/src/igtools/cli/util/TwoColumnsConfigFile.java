package igtools.cli.util;

import igtools.common.util.Pair;

import java.io.File;
import java.util.LinkedList;
import java.util.Map;
import java.util.Scanner;

public class TwoColumnsConfigFile {

	
	public static void parse(String file, Map<String,String> config) throws Exception{
		Scanner scanner = new Scanner(new File(file));
		String line;
		String[] cols;
		while(scanner.hasNextLine()){
			line = scanner.nextLine().trim();
			if(line.charAt(0) != '#'){
				cols = line.trim().split("\t");
				if(cols.length > 2 )
					throw new Exception("Too many columns, "+cols.length+", "+line);
				if(cols.length < 2)
					throw new Exception("Too less columns, "+cols.length+", "+line);
				config.put(cols[0], cols[1]);
			}
		}
		scanner.close();
	}
	
	
	public static void parse(String file, LinkedList<Pair<String,String>> config) throws Exception{
		Scanner scanner = new Scanner(new File(file));
		String line;
		String[] cols;
		while(scanner.hasNextLine()){
			line = scanner.nextLine().trim();
			if(line.charAt(0) != '#'){
				cols = line.trim().split("\t");
				if(cols.length > 2 )
					throw new Exception("Too many columns, "+cols.length+", "+line);
				if(cols.length < 2)
					throw new Exception("Too less columns, "+cols.length+", "+line);
				config.add(new Pair<String, String>(cols[0], cols[1]));
			}
		}
		scanner.close();
	}
}
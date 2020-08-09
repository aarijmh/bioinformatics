package com.course.bioinformatics;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Stream;

public class Week1 {
	
	public static HashSet<String> frequentPatterns(String text, int k)
	{
        		
		HashSet<String> frequentSets = new HashSet<>();
        List<Integer> countList = new ArrayList<>();
		
        int max = -1;
		for(int i =0 ; i < text.length() - k + 1 ; i++)
		{
			String pattern = text.substring(i,i+k);
			int count = countPattern(text,pattern);
			countList.add(count);
			if (count > max)
				max = count;
		}
		
		for(int i =0 ; i < text.length() - k + 1 ; i++)
		{
			if (countList.get(i) == max)
			{
				frequentSets.add(text.substring(i,i+k));
			}
		}
		
		return frequentSets;
	}

	public static int countPattern(String text, String pattern)
	{
		int count = 0;
		for (int i= 0; i< (text.length() - pattern.length() + 1);i++)
		{
			if(text.substring(i,i+pattern.length()).equals(pattern))
			{
				count++;
			}
		}
		return count;
	}
	
	public static List<Integer> indexPattern(String text, String pattern)
	{
		List<Integer> countList = new ArrayList<>();
		for (int i= 0; i< (text.length() - pattern.length() + 1);i++)
		{
			if(text.substring(i,i+pattern.length()).equals(pattern))
			{
				countList.add(i);
			}
		}
		return countList;
	}
	
	public static Set<String> clumpFinding2(String genome, int k, int l, int t)
	{
		long totalLength = genome.length();
		System.out.println("The length of the genome is : "+totalLength);
		Set<String> returnSet = new HashSet<>();
		Map<String,Integer> countMap = new HashMap<String, Integer>();
		
		for(int i = 0; i < l-k+1; i++)
		{
			String kmer = genome.substring(i,i+k);
			if(!countMap.containsKey(kmer))
			{
				countMap.put(kmer, 0);
			}
			int count = countMap.get(kmer);
			if(count + 1  >= t)
				returnSet.add(kmer);
			countMap.put(kmer, count+1);
		}
		
		for(int i = 1, j = l-k ; i < genome.length()- l - k+1 ;i++,j++ )
		{
				if(i % 100000 == 0)
					System.out.println(i + " Characters processed out of "+totalLength);
				
				String preKMER = genome.substring(i-1,i-1+k);
				countMap.put(preKMER, countMap.get(preKMER) - 1);
				
				String newKMER = genome.substring(j,j+k);
				if(!countMap.containsKey(newKMER))
				{
					countMap.put(newKMER, 0);
				}
				int count = countMap.get(newKMER);
				if(count + 1  >= t)
					returnSet.add(newKMER);
				countMap.put(newKMER, count+1);
		}
		
		System.out.println(returnSet.size());
		
		return returnSet;
	}
	
	public static Set<String> clumpFinding(String genome, int k, int l, int t)
	{
		Set<String> returnList = new HashSet<>();
		Map<String,Integer> countMap = new HashMap<String, Integer>();
		for(int i = 0; i < genome.length() - l; i++)
		{
			String subString = genome.substring(i,i+l);
			
			for (int j = 0 ; j < subString.length() - k + 1; j++)
			{
				String kmer = subString.substring(j,j+k);
				int c = countPattern(subString,kmer);
				if(c >= t)
				{
					if(!countMap.containsKey(kmer))
						countMap.put(kmer, 0);
					countMap.put(kmer, countMap.get(kmer) + 1);
				}
			}
			
			
		}
		return countMap.keySet();
	}
	static Map<String,String> complementHashMap = Map.of("A", "T", "T", "A", "G", "C", "C", "G","","");
	
	public static String reverseComplement(String text)
	{
		StringBuffer returnString = new StringBuffer();
		for(int i = text.length() - 1; i >= 0; i--)
		{
			returnString.append(complementHashMap.get(String.valueOf(text.charAt(i))));
		}
		return returnString.toString();
	}
	static Map<String,Integer> characterPosition = Map.of("A",0,"C",1,"G",2,"T",3);
	static Map<Integer,String> positionCharaceter = Map.of(0,"A",1,"C",2,"G",3,"T");
	
	public static Integer patternToInteger(String pattern)
	{
		int power = 0;
		Double sum  = 0.0;
		for(int j = pattern.length() - 1; j >= 0 ; j--)
		{
			sum += characterPosition.get(String.valueOf(pattern.charAt(j))) * Math.pow(4, power++);
			
		}
		return sum.intValue();
	}
	
	public static String numberToPattern(int position, int size)
	{
		List<String> pattern = new ArrayList<>();
		while(position > 0)
		{
			pattern.add(0,positionCharaceter.get(position % 4));
			position /= 4;
		}
		
		for(int i = pattern.size(); i < size;  i++)
		{
			pattern.add(0,"A");
		}
		
		return String.join("",pattern.subList(0, size));
	}
	
	public static List<Integer> computingFrequencies(String text, int k)
	{
		 List<Integer> frequencyList = new ArrayList<>();
		 
		 Double limit = Math.pow(4, k);
		 for(int i = 0; i < limit.intValue(); i++)
		 {
			 frequencyList.add(0);
		 }
		 
		 for(int i = 0; i < text.length() - k + 1; i++)
		 {
			 String pattern = text.substring(i,i+k);
			 int index = patternToInteger(pattern);
			 frequencyList.set(index, frequencyList.get(index) + 1);
		 }
		 
		 return frequencyList;
	}
	
	public static long patternToNumberRecursive(String pattern)
	{
		if (pattern.length() == 0)
			return 0;
		String symbol = String.valueOf(pattern.charAt(pattern.length()-1));
		pattern = pattern.substring(0,pattern.length()-1);
		return 4 * patternToNumberRecursive(pattern) +characterPosition.get(symbol);
	}
	

    
	public static void main(String[] args) {
		String text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT";
		String pattern  = "CTTGATCAT";
		int k = 3;
	
		
		System.out.println(countPattern("ACTGTACGATGATGTGTGTCAAAG", "TGT"));
		
		frequentPatterns(text,k).forEach(System.out::println);
		System.out.println(reverseComplement("TTGTGTC"));
		
//		System.out.println(reverseComplement(text));

		indexPattern("ATGACTTCGCTGTTACGCGC", "CGC").forEach(x->{
			System.out.print(x+" ");
		});
		
		/*
		 * StringBuilder contentBuilder = new StringBuilder();
		 * 
		 * try (Stream<String> stream = Files.lines(
		 * Paths.get("C:\\Users\\Aarij Mahmood\\Downloads\\Vibrio_cholerae.txt"),
		 * StandardCharsets.UTF_8)) { stream.forEach(s -> contentBuilder.append(s));
		 * 
		 * 
		 * System.out.println(contentBuilder.toString());
		 * 
		 * List<Integer> list = indexPattern(contentBuilder.toString(), pattern);
		 * 
		 * Path fileName = Path.of("e:\\demo.txt"); String content = "hello world !!";
		 * 
		 * StringBuffer buffer = new StringBuffer();
		 * list.forEach(x->buffer.append(x.toString()+" ")); Files.writeString(fileName,
		 * buffer.toString()); } catch (IOException e) { e.printStackTrace(); } try
		 * (Stream<String> stream = Files.lines(
		 * Paths.get("C:\\Users\\Aarij Mahmood\\Downloads\\E_coli.txt"),
		 * StandardCharsets.UTF_8)) { contentBuilder.setLength(0); stream.forEach(s ->
		 * contentBuilder.append(s));
		 * 
		 * String genome = contentBuilder.toString(); k = 9; int l = 500; int t = 3;
		 * long startTime = System.currentTimeMillis(); Set<String> kmers =
		 * clumpFinding2(genome, k, l, t); long endTime = System.currentTimeMillis();
		 * System.out.println("Total time (in ms) : "+(endTime - startTime));
		 * kmers.forEach(x->{ System.out.print(x + " "); }); } catch (Exception e) { //
		 * TODO Auto-generated catch block e.printStackTrace(); }
		 */
		
		/*
		 * System.out.println(patternToInteger("ATGCAA"));
		 * System.out.println(numberToPattern(7778, 7));
		 * 
		 * computingFrequencies(
		 * "CTATCCTGGATTATCTACGATTCATGTTCTGACACCCGCTAGCGAACAACTGATTAAGGCGTTGAAGAAGTAAGCTACACACAATTTACGGCCGACGCAACTTCGGGACCTCGTGAATTGCTAGTACTGCGAGATAATGTCAATGTGGCTGCCCTGGGGATATCACGTTCTGTGAGAATGCGAAATAGTCATGTGCTCGGCCCTAAGCCGGTATCGTGACCTAAAGGCGTGAATAGACCAAGTTGACAGGCGCTTGACATTAAAATGCAACAAGAGAGATATATCCCCGTACATTGGAAAGACAATGCATTGCCCATCGCTTGTTACCAAAAGTGGTTAACGAACCGCAAATGCAGCAGGTGTCGCCTGCCTGATTGCAAGCGTAGCTTGAGGGTGGCTTAATGGGCTCAGCTAATAGCTGTTCAGAACCTATATCAGGGGCCAAATTAGACTTACTCGGGTCCTATTCCGAGCGCGGGATATTAATAGGGTCAGTACGTACAAATTGAAGGACAGAGGTCTGGCAAGCCCAACAGGGCTATTACCCTGACCGACTCGGGTCACTCTTTGTGGGAGCAGGACCAATTAAGAGCTTATACAAGTCTTAGGAGGAAAGGTTTCATAGTAGAAATTTACACTGGTATAGGATTTCCATGACTCAGCTAGGTTAGCTCGTGCAATGTCACGCGCAGGACATTAACAGCGAAAACTCTTCGAGGGACGCCATCGGACCCGCCATTGTGATAGGGCAGTCCGTAGGGATGTCAAAAACGATCCTCTCAACT"
		 * ,6) .forEach(x->System.out.print(x+" ")); System.out.println();
		 * 
		 * System.out.println(patternToNumberRecursive("TGACTCTCAGGCTTTT"));
		 */
	}

}

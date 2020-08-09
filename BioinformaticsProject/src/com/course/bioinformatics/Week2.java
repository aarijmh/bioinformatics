package com.course.bioinformatics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

public class Week2 {
    
	public static List<Long> skewFinder(String genome)
	{
		List<Long> indices = new ArrayList<>();
		
		indices.add(0l);
		
		for(int i= 0; i < genome.length(); i++)
		{
			if(genome.charAt(i) == 'C')
			{
				indices.add(i + 1, indices.get(i)  - 1);
			}
			else if(genome.charAt(i) == 'G')
			{
				indices.add(i + 1, indices.get(i) + 1);
			}
			else {
				indices.add(i + 1, indices.get(i));
			}
		}
		
		return indices;
	}
	
	public static List<Long> skewMinimumFinder(String genome)
	{
		List<Long> indices = skewFinder(genome);
		
		List<Long> minIndices = new ArrayList<>();
		
		long minimum = Long.MAX_VALUE;
		
		long index = 0l;
		for(long l : indices)
		{
			if ( l < minimum)
			{
				minIndices.clear();
				minIndices.add(index);
				minimum = l;
			}
			else if( l == minimum)
			{
				minIndices.add(index);
			}
			index++;
		}
		
		return minIndices;
	}
	public static long hammingDistance(String str1, String str2)
	{
		long count = 0;
		for(int i = 0; i < str1.length() ; i ++)
		{
			if(str1.charAt(i) != str2.charAt(i))
				count++;
		}
		return count;
	}
	
	public static List<Long> approximateDistanceCalculator(String pattern, String text, int d)
	{
		List<Long> indices = new ArrayList<>();
		int patternLength = pattern.length();
		
		for(int i = 0; i < text.length() - pattern.length() + 1; i++)
		{
			long mismatches = hammingDistance(pattern, text.substring(i, i + patternLength));
			if(mismatches <= d)
			{
				indices.add((long)i);
			}
		}
		
		return indices;
	}
	
	public static Long approximateDistanceCalculatorCounter(String pattern, String text, int d)
	{
		Long count = 0l;
		int patternLength = pattern.length();
		
		for(int i = 0; i < text.length() - pattern.length() + 1; i++)
		{
			long mismatches = hammingDistance(pattern, text.substring(i, i + patternLength));
			if(mismatches <= d)
			{
				count++;
			}
		}
		
		return count;
	}
	
	public static HashSet<String> neighbours(String pattern, int d)
	{
		 HashSet<String> neighbourSet = new HashSet<>();
		 if(d  == 0)
		 {
			 neighbourSet.add(pattern);
		 }
		 else if(pattern.length() == 1)
		 {
			 neighbourSet.add("A");
			 neighbourSet.add("C");
			 neighbourSet.add("G");
			 neighbourSet.add("T");
		 }
		 else {
			 String suffixPattern = pattern.substring(1);
			 HashSet<String> suffixNeighbours = neighbours(suffixPattern,d);
			 for(String text : suffixNeighbours)
			 {
				 if(hammingDistance(suffixPattern, text) < d)
				 {
					 neighbourSet.add("A"+text);
					 neighbourSet.add("C"+text);
					 neighbourSet.add("G"+text);
					 neighbourSet.add("T"+text);
				 }
				 else {
					 neighbourSet.add(pattern.charAt(0)+text);
				 }
			 }
		 }
		 
		 return neighbourSet;
	}
	
	public static HashSet<String> frequentWordsWithMismatches(String text, int k, int d)
	{
		HashSet<String> frequentPatterns = new HashSet<>();
		TreeMap<Integer, Integer> sortedMap = new TreeMap<>();
		
		List<String> neighborhoods = new ArrayList<>();
		
		for (int i = 0; i < text.length() - k +1; i++)
		{
			neighborhoods.addAll(neighbours(text.substring(i,i+k), d));
		}
		
		int maximumValue = 1;
		for(int i = 0; i < neighborhoods.size(); i++)
		{
			String pattern = neighborhoods.get(i);
			int l =  Week1.patternToInteger(pattern);
			if(!sortedMap.containsKey(l))
			{
				sortedMap.put(l, 0);
			}
			int count = sortedMap.get(l) + 1;
			if(count > maximumValue)
				maximumValue = count;
			sortedMap.put(l,count);
		}
		
		
		for(Entry<Integer, Integer> entry : sortedMap.entrySet())
		{
			if(entry.getValue() == maximumValue)
			{
				frequentPatterns.add( Week1.numberToPattern(entry.getKey(), k));
			}
		}

		return frequentPatterns;
	}
	
	public static HashSet<String> frequentWordsWithMismatchesWithReverseComplement(String text, int k, int d)
	{
		HashSet<String> frequentPatterns = new HashSet<>();
		TreeMap<Integer, Integer> sortedMap = new TreeMap<>();
		
		List<String> neighborhoods = new ArrayList<>();
		
		for (int i = 0; i < text.length() - k +1; i++)
		{
			String pattern = text.substring(i,i+k);
			String complement = Week1.reverseComplement(pattern);
			neighborhoods.addAll(neighbours(pattern, d));
			neighborhoods.addAll(neighbours(complement, d));
		}
		
		int maximumValue = 1;
		for(int i = 0; i < neighborhoods.size(); i++)
		{
			String pattern = neighborhoods.get(i);
			int l =  Week1.patternToInteger(pattern);
			if(!sortedMap.containsKey(l))
			{
				sortedMap.put(l, 0);
			}
			int count = sortedMap.get(l) + 1;
			if(count > maximumValue)
				maximumValue = count;
			sortedMap.put(l,count);
		}
		
		
		for(Entry<Integer, Integer> entry : sortedMap.entrySet())
		{
			if(entry.getValue() == maximumValue)
			{
				frequentPatterns.add( Week1.numberToPattern(entry.getKey(), k));
			}
		}

		return frequentPatterns;
	}

	
	public static void main(String[] args) {
		//1.3.1
		skewFinder("CATTCCAGTACTTCATGATGGCGTGAAGA").forEach(x->System.out.print(x+ " "));
		System.out.println();
		//1.3.2
//		skewMinimumFinder("CATTCCAGTACTTCATGATGGCGTGAAGA").forEach(x->System.out.print(x+ " "));

		//1.3.3
//		System.out.println(hammingDistance("CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT",
//				                           "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"));
		
		//1.3.4
//		List<String> stringList = UtilityFunction.readStringListFromFile("C:\\Users\\Aarij Mahmood\\Downloads\\dataset_9_4.txt");
//		approximateDistanceCalculator("CCC",
//									  "CATGCCATTCGCATTGTCCCAGTGA",
//									  2)
//									.forEach(x->System.out.print(x +" "));
		
		//1.3.5
		System.out.println(approximateDistanceCalculatorCounter("CCC",
									  "CATGCCATTCGCATTGTCCCAGTGA",
									  2));
		
		//Test neighbours
//		neighbours("ACGT", 4).forEach(x->System.out.print(x + " "));
//		System.out.println(neighbours("ACGT", 4).size());
		//Test frequentWordsWithMismatches
		frequentWordsWithMismatches("ACGT", 4,3).forEach(x->System.out.print(x + " "));
		System.out.println();
		System.out.println(frequentWordsWithMismatches("ACGT", 4,3).size());
		
		
		//Test frequentWordsWithMismatches with complement
//		frequentWordsWithMismatchesWithReverseComplement("TGTTGTCTCCTCGATGTGCAGAGCAGGTGAGGTGCAGGTCTCGGTGGTGAGATGTTGTGGTGCAGAGCAGACTCTGTGGTGAGAGATGTGAGAGGTGACTCGGTTGTGCATGTCTCCTCGAGGTGGTTGTCTCCTCGCAGGTGATGTCTCGGTGCAGAGCAGGTGGTGAGCACTCGCAGGTCTCCTCGCAGGTCTCGGTTGTGGTGGTGCACTCGAGCAGCACTCGCATGTGCA", 6,3).forEach(x->System.out.print(x + " "));
	}

}

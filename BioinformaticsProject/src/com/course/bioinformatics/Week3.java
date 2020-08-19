package com.course.bioinformatics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


public class Week3 {
	

	public static ArrayList<String>  motifEnumeration(List<String> dna, Integer k, Integer d)
	{
		java.util.HashSet<String> patterns = new java.util.HashSet<>();
		
		String firstPattern = dna.get(0);
		
		for(int i = 0; i < firstPattern.length() - k + 1; i++)
		{
			String subPattern = firstPattern.substring(i,i+k);
			java.util.HashSet<String> neighbours = Week2.neighbours(subPattern, d);
			
			for(String subPat : neighbours)
			{
				boolean match = false;
				for(int j = 1; j < dna.size(); j++)
				{
					match = false;
					String secondPattern = dna.get(j);
					for(int l = 0; l < secondPattern.length() - k + 1; l++)
					{
						if(Week2.hammingDistance(subPat, secondPattern.substring(l,l+k) ) <= d)
						{
							match = true;
							break;
						}
					}
					if (!match)
						break;
				}
				if(match)
					patterns.add(subPat);
				else
					continue;
			}
		}
		
		return new ArrayList<String>(patterns);
	}

	public static int distanceBetweenPatternAndStrings(String pattern, List<String> dna)
	{
		int k = pattern.length();
		int distance = 0;
		
		for(String text : dna)
		{
			int hammingDistance = Integer.MAX_VALUE;
			
			for(int i = 0; i < text.length() - k + 1; i++)
			{
				String patternPrime = text.substring(i,i+k);
				int temp = (int) Week2.hammingDistance(pattern, patternPrime);
				if(hammingDistance > temp)
				{
					hammingDistance = temp;
				}
			}
			distance += hammingDistance;
		}
		
		return distance;
	}
	
	public static String medianString(List<String> dna, int k)
	{
		String median = "";
		int distance = Integer.MAX_VALUE;
		for (int i = 0; i < Math.pow(4, k);i++)
		{
			String pattern = Week1.numberToPattern(i, k);
			int temp = distanceBetweenPatternAndStrings(pattern, dna);
			
			if(temp == 0)
				System.out.println(pattern);
			if( distance > temp)
			{
				distance = temp;
				median = pattern;
			}
		}
		System.out.println(distance);
		return median;
	}
	public static String profileMostProbablekmer(String text, int k, Map<Character, List <Double>> profile) {
		
		double distance = -1.0;
		String median = "";
		for(int i = 0; i < text.length() - k +1 ;i++)
		{
			double tempDistance = 1.0;
			String pattern = text.substring(i,i+k);
			for(int j = 0; j < pattern.length(); j++)
			{
				Character c = pattern.charAt(j);
				tempDistance *= profile.get(c).get(j);
			}
			if(tempDistance > distance)
			{
				distance = tempDistance;
				median = pattern;
			}
		}

		return median;
	}
	public static List<Double> createEmptyList(int m)
	{
		List<Double> dd = new ArrayList<>();
		for(int i=0;i<m;i++)
			dd.add(0.0);
		return dd;
	}
	
	public static int scoreMotifs(List<String> motifs,int k)
	{
		int score = 0;
		for(int i = 0 ; i < k ; i++)
		{
			java.util.Map<Character, Integer> positionMap = new java.util.HashMap<>();
			positionMap.put('A', 0);
			positionMap.put('C', 0);
			positionMap.put('G', 0);
			positionMap.put('T', 0);
			
			for(int j  = 0; j < motifs.size();j++)
			{
				positionMap.put(motifs.get(j).charAt(i), positionMap.get(motifs.get(j).charAt(i)) + 1);
			}
			
			int maxFrequency = -1;
			
			for(Integer c : positionMap.values())
			{
				if(c > maxFrequency)
					maxFrequency = c;
			}
			score += motifs.size() - maxFrequency;
		}
		
		return score;
	}
	public static List <String> greedyMotifSearch(List <String> dna, int k, int t) {
		java.util.Map<Character, Integer> positionMap = new java.util.HashMap<>();
		positionMap.put('A', 0);
		positionMap.put('C', 1);
		positionMap.put('G', 2);
		positionMap.put('T', 3);
		
		
		List<String> bestMotifs = new ArrayList<>();
		for(int i = 0; i < dna.size(); i++)
		{
			bestMotifs.add(dna.get(i).substring(0,k));
		}
		
		String firstString = dna.get(0);
		for(int i = 0; i < firstString.length() - k + 1; i++)
		{
			String motif = firstString.substring(i,i+k);
			List<String> motifList = new ArrayList<>();
			motifList.add(motif);
			for(int j = 1; j < t; j++)
			{
				java.util.Map<Character,List<Double>> profile = new java.util.HashMap<>();
				profile.put('A', createEmptyList(k));
				profile.put('C', createEmptyList(k));
				profile.put('G', createEmptyList(k));
				profile.put('T', createEmptyList(k));

				for(int m = 0; m < k; m++)
				{
					
					java.util.Map<Character,Integer> frequencyMap = new java.util.LinkedHashMap<>();
					frequencyMap.put('A', 0);
					frequencyMap.put('C', 0);
					frequencyMap.put('G', 0);
					frequencyMap.put('T', 0);
					for(int n = 0 ; n < motifList.size() ; n++)
					{
						frequencyMap.put(motifList.get(n).charAt(m), frequencyMap.get(motifList.get(n).charAt(m)) + 1);
					}
					profile.get('A').set(m, frequencyMap.get('A')/(double)motifList.size());
					profile.get('C').set(m, frequencyMap.get('C')/(double)motifList.size());
					profile.get('G').set(m, frequencyMap.get('G')/(double)motifList.size());
					profile.get('T').set(m, frequencyMap.get('T')/(double)motifList.size());
				}
				motifList.add(profileMostProbablekmer(dna.get(j), k,profile));
			}
			if(scoreMotifs(motifList, k) < scoreMotifs(bestMotifs, k))
			{
				bestMotifs = new ArrayList<>(motifList);
			}
		}
		
		return bestMotifs;
	}
	
	
	public static List <String> greedyMotifSearchPseudoCounts(List <String> dna, int k, int t) {
		java.util.Map<Character, Integer> positionMap = new java.util.HashMap<>();
		positionMap.put('A', 0);
		positionMap.put('C', 1);
		positionMap.put('G', 2);
		positionMap.put('T', 3);
		
		
		List<String> bestMotifs = new ArrayList<>();
		for(int i = 0; i < dna.size(); i++)
		{
			bestMotifs.add(dna.get(i).substring(0,k));
		}
		
		String firstString = dna.get(0);
		for(int i = 0; i < firstString.length() - k + 1; i++)
		{
			String motif = firstString.substring(i,i+k);
			List<String> motifList = new ArrayList<>();
			motifList.add(motif);
			for(int j = 1; j < t; j++)
			{
				java.util.Map<Character,List<Double>> profile = new java.util.HashMap<>();
				profile.put('A', createEmptyList(k));
				profile.put('C', createEmptyList(k));
				profile.put('G', createEmptyList(k));
				profile.put('T', createEmptyList(k));

				for(int m = 0; m < k; m++)
				{
					
					java.util.Map<Character,Integer> frequencyMap = new java.util.LinkedHashMap<>();
					frequencyMap.put('A', 1);
					frequencyMap.put('C', 1);
					frequencyMap.put('G', 1);
					frequencyMap.put('T', 1);
					for(int n = 0 ; n < motifList.size() ; n++)
					{
						frequencyMap.put(motifList.get(n).charAt(m), frequencyMap.get(motifList.get(n).charAt(m)) + 1);
					}
					profile.get('A').set(m, frequencyMap.get('A')/((double)motifList.size()+4.0));
					profile.get('C').set(m, frequencyMap.get('C')/((double)motifList.size()+4.0));
					profile.get('G').set(m, frequencyMap.get('G')/((double)motifList.size()+4.0));
					profile.get('T').set(m, frequencyMap.get('T')/((double)motifList.size()+4.0));
				}
				motifList.add(profileMostProbablekmer(dna.get(j), k,profile));
			}
			if(scoreMotifs(motifList, k) < scoreMotifs(bestMotifs, k))
			{
				bestMotifs = new ArrayList<>(motifList);
			}
		}
		
		return bestMotifs;
	}

    		
	public static void main(String[] args) {
		List<String> input = UtilityFunction.readStringListFromFile("input.txt");
		
		List<Character> characters = List.of('A','C','G','T');
		
		
		// Test motifEnumeration
//		String [] inp = input.get(0).split(" ");
//		
//		motifEnumeration(input.subList(1,input.size()),Integer.valueOf(inp[0]),Integer.valueOf(inp[1]))
//		.forEach(x->System.out.print(x+" "));
		
		//Test distanceBetweenPatternAndStrings
		
//		
//		System.out.println(distanceBetweenPatternAndStrings(input.get(0),List.of(input.get(1).split(" "))));
		
		//Test medianString
		
//		System.out.println(medianString(input.subList(1, input.size()), Integer.valueOf(input.get(0))));
		
		//Test profileMostProbablekmer
//		Map<Character,List<Double>> profile = new HashMap<>();
//		String text = input.get(0);
//		int k = Integer.valueOf(input.get(1));
//		for(int i = 0; i < characters.size(); i++)
//		{
//			String [] st = input.get(i + 2).split(" ");
//			List<Double> pl = new ArrayList<>();
//			for(int j = 0; j < st.length; j++)
//			{
//				pl.add(Double.valueOf(st[j]));
//			}
//			profile.put(characters.get(i), pl);
//		}
//		System.out.println(profileMostProbablekmer(text,k,profile));
		
		
		//Test greedy motif search
//		String [] inputFirst = input.get(0).split(" ");
//		greedyMotifSearch(input.subList(1, input.size()), Integer.valueOf(inputFirst[0]), Integer.valueOf(inputFirst[1]))
//		.forEach(x->System.out.print(x + " "));
		
//		Test greedy motif search pseudocount
//		String [] inputFirst = input.get(0).split(" ");
//		greedyMotifSearchPseudoCounts(input.subList(1, input.size()), Integer.valueOf(inputFirst[0]), Integer.valueOf(inputFirst[1]))
//		.forEach(x->System.out.println(x));
		
		List<String> testList = List.of("CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
				"GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
				"GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG");
		System.out.println(medianString(testList, 7));
		
//		
		Map<Character, List<Double>> profile = 
		 Map.of('A',List.of(0.4, 0.3, 0.0, 0.1, 0.0, 0.9),
				'C',List.of(0.2, 0.3, 0.0, 0.4, 0.0, 0.1),
				'G',List.of(0.1, 0.3, 1.0, 0.1, 0.5, 0.0),
				'T',List.of(0.3, 0.1, 0.0, 0.4, 0.5, 0.0));
		
		double tempDistance = 1.0;
		String pattern ="CAGTGA";
		for(int j = 0; j < pattern.length(); j++)
		{
			Character c = pattern.charAt(j);
			tempDistance *= profile.get(c).get(j);
		}
		System.out.println("Profile  : " +  tempDistance);
	}

}

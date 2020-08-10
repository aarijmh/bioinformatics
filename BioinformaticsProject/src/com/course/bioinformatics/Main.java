package com.course.bioinformatics;

import java.util.ArrayList;
import java.util.List;

public class Main {
	//fill in your GreedyMotifSearchWithPseudocounts() function here with any subroutines you need.

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
	public static String profileMostProbablekmer(String text, int k, java.util.Map<Character, List <Double>> profile) {
			
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

	public static List <String> GreedyMotifSearchWithPseudocounts(List <String> dna, int k, int t) {
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

}

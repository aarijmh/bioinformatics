package com.course.bioinformatics;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Week4 {
	static Random random = new Random();

	public static List<Integer> generateKRandomIntsinRange(int k, int size) {
		List<Integer> randomInts = new ArrayList<Integer>();
		for (int i = 0; i < k; i++)
			randomInts.add(random.nextInt(size));
		return randomInts;
	}

	public static List<Double> createEmptyList(int m) {
		List<Double> dd = new ArrayList<>();
		for (int i = 0; i < m; i++)
			dd.add(0.0);
		return dd;
	}

	public static java.util.Map<Character, List<Double>> createProfile(List<String> motifList, int k) {
		java.util.Map<Character, List<Double>> profile = new java.util.HashMap<>();
		profile.put('A', createEmptyList(k));
		profile.put('C', createEmptyList(k));
		profile.put('G', createEmptyList(k));
		profile.put('T', createEmptyList(k));

		for (int m = 0; m < k; m++) {

			java.util.Map<Character, Integer> frequencyMap = new java.util.LinkedHashMap<>();
			frequencyMap.put('A', 1);
			frequencyMap.put('C', 1);
			frequencyMap.put('G', 1);
			frequencyMap.put('T', 1);
			for (int n = 0; n < motifList.size(); n++) {
				frequencyMap.put(motifList.get(n).charAt(m), frequencyMap.get(motifList.get(n).charAt(m)) + 1);
			}
			profile.get('A').set(m, frequencyMap.get('A') / ((double) motifList.size() + 4.0));
			profile.get('C').set(m, frequencyMap.get('C') / ((double) motifList.size() + 4.0));
			profile.get('G').set(m, frequencyMap.get('G') / ((double) motifList.size() + 4.0));
			profile.get('T').set(m, frequencyMap.get('T') / ((double) motifList.size() + 4.0));
		}

		return profile;
	}

	public static List<String> getMotifs(List<String> dna, java.util.Map<Character, List<Double>> profile, int k) {
		List<String> motifs = new ArrayList<>();
		for (int i = 0; i < dna.size(); i++) {
			motifs.add(Week3.profileMostProbablekmer(dna.get(i), k, profile));
		}
		return motifs;
	}

	public static List<String> randomizedMotifSearch(List<String> dna, int k, int t) {
		List<String> bestMotifs = new ArrayList<>();

//		for (int i = 0; i < dna.size(); i++) {
//			int index = random.nextInt(dna.get(i).length() - k);
//		bestMotifs.add(dna.get(i).substring(index, index + k));
//		}
		bestMotifs.add("CCA");
		bestMotifs.add("CCT");
		bestMotifs.add("CTT");
		bestMotifs.add("TTG");
		
		List<String> motifList = new ArrayList<>();
		while (true) {
			java.util.Map<Character, List<Double>> profile = createProfile(bestMotifs, k);
			motifList = getMotifs(dna, profile, k);
			if (Week3.scoreMotifs(motifList, k) < Week3.scoreMotifs(bestMotifs, k)) {
				bestMotifs = motifList;
			} else {
				return motifList;
			}
		}

	}
	
	private static List<String> removeKthIndexFromList(List<String> l, int index)
	{
		 List<String> returnList = new ArrayList<String>();
		 for(String s : l)
			 returnList.add(s);
		 returnList.remove(index);
		 return returnList;
	}
	
	public static int mostProbableKIndex(List<Double> probabilities, int n)
	{
		double sum = 0;
		for(Double d : probabilities)
			sum += d;
		
		double a = random.nextDouble() * sum;
		double sum_pr = 0.0;
		for(int  i =0; i < n-1 ;i++)
		{
			if (a >= sum_pr && a < (sum_pr + probabilities.get(i)))
			{
				return i;
			}
			sum_pr += probabilities.get(i);
		}
		return n-1;
	}

	
	public static String calculateWeightedKMer(String text, java.util.Map<Character, List<Double>> profile, int k)
	{
		List<String> possibleKMers = new ArrayList<>();
		for(int i = 0; i< text.length() - k +1 ; i++)
		{
			possibleKMers.add(text.substring(i,i+k));
		}
		
		List<Double> probabilities = new ArrayList<>();
		
		for(String kmer : possibleKMers)
		{
			double prob = 1.0;
			for(int j = 0; j < kmer.length(); j++)
				prob *= profile.get(kmer.charAt(j)).get(j);
			probabilities.add(prob);
		}
		
		return possibleKMers.get(mostProbableKIndex(probabilities,text.length() - k + 1));
	}

	public static List<String> gibbsRandomizedMotifSearch(List<String> dna, int k, int t, int n) {
		List<String> bestMotifs = new ArrayList<>();

		for (int i = 0; i < dna.size(); i++) {
			int index = random.nextInt(dna.get(i).length() - k);
			bestMotifs.add(dna.get(i).substring(index, index + k));
		}

		for (int i = 0; i < n; i++) {
			
			int index = random.nextInt(t);
			List<String> motifList = new ArrayList<>();

				java.util.Map<Character, List<Double>> profile = createProfile(removeKthIndexFromList(bestMotifs,index), k);
				String motif = calculateWeightedKMer(dna.get(index),profile,k);
				
				motifList = getMotifs(dna, profile, k);
				motifList.set(index, motif);
				if (Week3.scoreMotifs(motifList, k) < Week3.scoreMotifs(bestMotifs, k)) {
					bestMotifs = motifList;
				} 

		}
		return bestMotifs;
	}



	public static void main(String[] args) {
		List<String> input = UtilityFunction.readStringListFromFile("input.txt");

		List<Character> characters = List.of('A', 'C', 'G', 'T');

		// Test randomizedMotifSearch

//		String [] i = input.get(0).split(" ");
//		
//		List<String> lastMotifs;
//		int count = 0;
//		int k = Integer.valueOf(i[0]);
//		int t = Integer.valueOf(i[1]);
//		lastMotifs = randomizedMotifSearch(input.subList(1, input.size()),k,t);
//
//		while(count < 1000)
//		{
//			List<String> bestmotifs = randomizedMotifSearch(input.subList(1, input.size()),k,t);
//			if(Week3.scoreMotifs(bestmotifs, k) < Week3.scoreMotifs(lastMotifs, k))
//				lastMotifs = bestmotifs;
//			count++;
//		}
//		
//		lastMotifs.forEach(System.out::println);

//		double p1 = (600.0 - 15.0) / (600.0 - 15.0 + 1.0);
//		double pC = Math.pow(p1, 10);
//		double pAnswer = 1 - pC;
//		System.out.println(pAnswer);
		
//		String [] i = input.get(0).split(" ");
//		
//		String [] i = input.get(0).split(" ");
//		List<String> lastMotifs;
//		int count = 0;
//		int k = Integer.valueOf(i[0]);
//		int t = Integer.valueOf(i[1]);
//		int n =  Integer.valueOf(i[2]);
		
//				lastMotifs = gibbsRandomizedMotifSearch(input.subList(1, input.size()),k,t,n);
//	
//				while(count < 20)
//				{
//					long start = System. currentTimeMillis();
//					List<String> bestmotifs = gibbsRandomizedMotifSearch(input.subList(1, input.size()),k,t,n);
//					if(Week3.scoreMotifs(bestmotifs, k) < Week3.scoreMotifs(lastMotifs, k))
//						lastMotifs = bestmotifs;
//					long end = System. currentTimeMillis();
//					count++;
//					System.out.println(count +  "  ---  "+(end-start));
//				}
//				
//				lastMotifs.forEach(System.out::println);
		
		List<String> inputList = new ArrayList<>();
		inputList.add("AAGCCAAA");
		inputList.add("AATCCTGG");
		inputList.add("GCTACTTG");
		inputList.add("ATGTTTTG");
		randomizedMotifSearch(inputList,3,4).forEach(x->System.out.print(x+","));;

	}
}

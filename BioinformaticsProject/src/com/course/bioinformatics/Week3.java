package com.course.bioinformatics;

import java.util.ArrayList;
import java.util.List;

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
			if( distance > temp)
			{
				distance = temp;
				median = pattern;
			}
		}
		return median;
	}

    		
	public static void main(String[] args) {
		// Test motifEnumeration
//		List<String> input = UtilityFunction.readStringListFromFile("input.txt");
//		String [] inp = input.get(0).split(" ");
//		
//		motifEnumeration(input.subList(1,input.size()),Integer.valueOf(inp[0]),Integer.valueOf(inp[1]))
//		.forEach(x->System.out.print(x+" "));
		
		//Test distanceBetweenPatternAndStrings
		
//		List<String> input = UtilityFunction.readStringListFromFile("input.txt");
//		
//		System.out.println(distanceBetweenPatternAndStrings(input.get(0),List.of(input.get(1).split(" "))));
		
		//Test medianString
		
		List<String> input = UtilityFunction.readStringListFromFile("input.txt");
		
		System.out.println(medianString(input.subList(1, input.size()), Integer.valueOf(input.get(0))));

	}

}

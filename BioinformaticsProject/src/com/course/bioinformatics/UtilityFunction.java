package com.course.bioinformatics;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class UtilityFunction {

	public static String readFromFile(String path) {
		StringBuffer contentBuilder = new StringBuffer();
		try (Stream<String> stream = Files.lines(Paths.get(path), StandardCharsets.UTF_8)) {
			stream.forEach(s ->{ 
				contentBuilder.append(s);
				});

		} catch (Exception e) {
			e.printStackTrace();
		}
		return contentBuilder.toString();
	}
	
	public static List<String> readStringListFromFile(String path) {
		List<String> returnList = new ArrayList<>();
		try (Stream<String> stream = Files.lines(Paths.get(path), StandardCharsets.UTF_8)) {
			stream.forEach(s ->{ 
				returnList.add(s);
				});

		} catch (Exception e) {
			e.printStackTrace();
		}
		return returnList;
	}
	
	public static void quickSort(int arr[], int begin, int end) {
	    if (begin < end) {
	        int partitionIndex = partition(arr, begin, end);
	 
	        quickSort(arr, begin, partitionIndex-1);
	        quickSort(arr, partitionIndex+1, end);
	    }
	}
	
	private static int partition(int arr[], int begin, int end) {
	    int pivot = arr[end];
	    int i = (begin-1);
	 
	    for (int j = begin; j < end; j++) {
	        if (arr[j] <= pivot) {
	            i++;
	 
	            int swapTemp = arr[i];
	            arr[i] = arr[j];
	            arr[j] = swapTemp;
	        }
	    }
	 
	    int swapTemp = arr[i+1];
	    arr[i+1] = arr[end];
	    arr[end] = swapTemp;
	 
	    return i+1;
	}
}

/** 
*Used to outsource the adjacency matrix calculations from MatLab
*@Author George Emanuel
*emanuega0@gmail.com
*Copyright Presidents and Fellows of Harvard College, 2017.
*/
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


public class AdjacencyMatrix {
	
	
	public static void main(String[] args) {
		int[][] test = {{1, 0}, {0, 1}};
		calculateEdges(test, 1);
	}
	
	public static int[][] calculateEdges(int[][] input, int cutOff) {
		
		ArrayList<LinkedList<Integer>> edges = new ArrayList<LinkedList<Integer>>(input.length);
		
		for (int i=0;i<input.length;i++) {
			LinkedList<Integer> currentEdges = new LinkedList<Integer>();
			for (int j=i+1;j<input.length;j++) {
				int currentDistance = distance(input[i], input[j], cutOff);
				if (currentDistance <= cutOff)
					currentEdges.add(new Integer(j));
			}
			edges.add(i, currentEdges);
			if (i % 1000 == 0)
				System.out.println(i);
			
		}
		
		int edgeCount = 0;
		for (int i=0;i<input.length;i++) {
			edgeCount += edges.get(i).size();
		}
		
		int[][] edgeMatrix = new int[edgeCount][2];
		int edgeIndex = 0;
		for (int i=0;i<edges.size();i++) {
			for (Integer currentTo : edges.get(i)) {
				edgeMatrix[edgeIndex][0] = i;
				edgeMatrix[edgeIndex][1] = currentTo.intValue();
				edgeIndex++;
			}
		}
		
		
		return edgeMatrix;
	}
	
	/** Find the nearest element in sequenceLibrary that is closest to each of the elements in 
	* inputSequence
	*/
	public static int[] calculateNearestSequences(int[][] inputSequences, int[][] sequenceLibrary) {
		
		int[] nearest = new int[inputSequences.length];
		
		for (int i=0;i<inputSequences.length;i++) {
			nearest[i] = nearestSequenceFromLibrary(inputSequences[i], sequenceLibrary);
		}
		
		return nearest;
	}
	
	private static int nearestSequenceFromLibrary(int[] inputSequence, int[][] sequenceLibrary) {
		int minDistance = inputSequence.length+1;
		int minIndex = -1;
		
		for (int i=0;i<sequenceLibrary.length;i++) {
			int currentDistance = distance(inputSequence, sequenceLibrary[i], minDistance);
			if (currentDistance < minDistance) {
				minDistance = currentDistance;
				minIndex = i;
			}
		}
		
		return minIndex;
	}
	
	
	private static int distance(int[] input1, int[] input2, int cutOff) {
		int distance = 0;
		for (int i=0;i<input1.length;i++) {
			if (input1[i] != input2[i]) {
				distance++;
				if (distance > cutOff) return distance;
			}
		}
		return distance;
	}
}

/***********************************************************************
* FILENAME :        COMP20003_668329.c
*
* DESCRIPTION :
*       Dijkstra shortest path algorithm
*
* SUBJECT :
*       COMP20003 
* 
* AUTHOR :    	Rauan Kelesbekov        
* START DATE :	11 Apr 2014
* CONTACTS :  	rkelesbekov@student.unimelb.edu.au
* STUDENT ID : 	668329
* CHANGES : 	Code was modified for 
* 				Operations Research subject on 12 Nov 2020
* 				Removed Malloc and added main function
*
************************************************************************/

#include <stdio.h>
#define INFINITY 9999999
#define maxVertices 200

void Dijkstra(int adjMX[maxVertices][maxVertices], int n, int initial, int end);
void Dijkstra(int adjMX[maxVertices][maxVertices], int n, int initial, int end) {
	int cost[maxVertices][maxVertices], distance[maxVertices], preceding[maxVertices];
	int visited[maxVertices], count, minDist, lastNode, i, j;
	/* Fill your function algorithm below.
		malloc and output requirements are listed in your assignment sheet.
	*/

	// cost matrix gen
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (adjMX[i][j] == 0){
				cost[i][j] = INFINITY; // fill with heavy price
			}
			else{
				cost[i][j] = adjMX[i][j]; // use the adjMX transport cost 
			}
		}
	}
	for (i = 0; i < n; i++) {
		distance[i] 	= cost[initial][i]; // distance list
		preceding[i] 	= initial;			// preceding list
		visited[i] 		= 0; 				// visited list
	}

	distance[initial] 	= 0;	// least dist to vert
	visited[initial] 	= 1; 	// 0 for unvisited, 1 for visited
	count 				= 1;	// nu of vertices visited

	while (count < n - 1) { // untill all vertices visited
		minDist = INFINITY;
		for (i = 0; i < n; i++){
			if (distance[i] < minDist && !visited[i]) {
				minDist = distance[i];
				lastNode = i;
			}
		}

		visited[lastNode] = 1;
		for (i = 0; i < n; i++){
			if (!visited[i]){
				if (minDist + cost[lastNode][i] < distance[i]) {
					distance[i] = minDist + cost[lastNode][i]; // total cost to get to the vertex form start
					preceding[i] = lastNode;	// previous vertex to parse the path
				}
			}
		}
		count++;
	}

	
	// out
	printf("  to | cost |  via"); 
		for (i = 0; i < n; i++){
			if (i != initial) {
				printf("\n%4d | %4d | %4d", i, distance[i], preceding[i]);
			}
		}
		printf("\nWhere:\n\"to\"\t- the destination vertex");
		printf("\n\"cost\"\t- total cost of transportation to the vertex from initial");
		printf("\n\"via\"\t- preceding vertex to find the path\n");
		int last = preceding[end];
		printf("Path: end ");
		while (last != initial){
			printf("%3d <-", last);
			last = preceding[last];
		}
		printf(" start\n");
	}
	
// leave main field empty

// TEST FUNCTION
int main() {
	int adjMX[maxVertices][maxVertices], i, j, n, initial, end;
	n = 12;

	adjMX[0][0] = 0;
	adjMX[0][1] = 20;
	adjMX[0][2] = 15;
	adjMX[0][3] = 0;
	adjMX[0][4] = 0;
	adjMX[0][5] = 0;
	adjMX[0][6] = 0;
	adjMX[0][7] = 0;
	adjMX[0][8] = 0;
	adjMX[0][9] = 0;
	adjMX[0][10] = 0;
	adjMX[0][11] = 0;

	adjMX[1][0] = 0;
	adjMX[1][1] = 0;
	adjMX[1][2] = 0;
	adjMX[1][3] = 12;
	adjMX[1][4] = 32;
	adjMX[1][5] = 0;
	adjMX[1][6] = 0;
	adjMX[1][7] = 0;
	adjMX[1][8] = 0;
	adjMX[1][9] = 0;
	adjMX[1][10] = 0;
	adjMX[1][11] = 0;

	adjMX[2][0] = 0;
	adjMX[2][1] = 0;
	adjMX[2][2] = 0;
	adjMX[2][3] = 28;
	adjMX[2][4] = 0;
	adjMX[2][5] = 0;
	adjMX[2][6] = 35;
	adjMX[2][7] = 0;
	adjMX[2][8] = 0;
	adjMX[2][9] = 0;
	adjMX[2][10] = 0;
	adjMX[2][11] = 0;

	adjMX[3][0] = 0;
	adjMX[3][1] = 0;
	adjMX[3][2] = 0;
	adjMX[3][3] = 0;
	adjMX[3][4] = 17;
	adjMX[3][5] = 45;
	adjMX[3][6] = 0;
	adjMX[3][7] = 45;
	adjMX[3][8] = 0;
	adjMX[3][9] = 0;
	adjMX[3][10] = 0;
	adjMX[3][11] = 0;

	adjMX[4][0] = 0;
	adjMX[4][1] = 0;
	adjMX[4][2] = 0;
	adjMX[4][3] = 0;
	adjMX[4][4] = 0;
	adjMX[4][5] = 12;
	adjMX[4][6] = 0;
	adjMX[4][7] = 0;
	adjMX[4][8] = 0;
	adjMX[4][9] = 0;
	adjMX[4][10] = 0;
	adjMX[4][11] = 0;

	adjMX[5][0] = 0;
	adjMX[5][1] = 0;
	adjMX[5][2] = 0;
	adjMX[5][3] = 0;
	adjMX[5][4] = 0;
	adjMX[5][5] = 0;
	adjMX[5][6] = 0;
	adjMX[5][7] = 11;
	adjMX[5][8] = 17;
	adjMX[5][9] = 0;
	adjMX[5][10] = 0;
	adjMX[5][11] = 0;

	adjMX[6][0] = 0;
	adjMX[6][1] = 0;
	adjMX[6][2] = 0;
	adjMX[6][3] = 0;
	adjMX[6][4] = 0;
	adjMX[6][5] = 0;
	adjMX[6][6] = 0;
	adjMX[6][7] = 9;
	adjMX[6][8] = 0;
	adjMX[6][9] = 12;
	adjMX[6][10] = 0;
	adjMX[6][11] = 0;

	adjMX[7][0] = 0;
	adjMX[7][1] = 0;
	adjMX[7][2] = 0;
	adjMX[7][3] = 0;
	adjMX[7][4] = 0;
	adjMX[7][5] = 0;
	adjMX[7][6] = 0;
	adjMX[7][7] = 0;
	adjMX[7][8] = 8;
	adjMX[7][9] = 0;
	adjMX[7][10] = 5;
	adjMX[7][11] = 0;

	adjMX[8][0] = 0;
	adjMX[8][1] = 0;
	adjMX[8][2] = 0;
	adjMX[8][3] = 0;
	adjMX[8][4] = 0;
	adjMX[8][5] = 0;
	adjMX[8][6] = 0;
	adjMX[8][7] = 0;
	adjMX[8][8] = 0;
	adjMX[8][9] = 0;
	adjMX[8][10] = 0;
	adjMX[8][11] = 7;

	adjMX[9][0] = 0;
	adjMX[9][1] = 0;
	adjMX[9][2] = 0;
	adjMX[9][3] = 0;
	adjMX[9][4] = 0;
	adjMX[9][5] = 0;
	adjMX[9][6] = 0;
	adjMX[9][7] = 0;
	adjMX[9][8] = 0;
	adjMX[9][9] = 0;
	adjMX[9][10] = 0;
	adjMX[9][11] = 18;

	adjMX[10][0] = 0;
	adjMX[10][1] = 0;
	adjMX[10][2] = 0;
	adjMX[10][3] = 0;
	adjMX[10][4] = 0;
	adjMX[10][5] = 0;
	adjMX[10][6] = 0;
	adjMX[10][7] = 0;
	adjMX[10][8] = 0;
	adjMX[10][9] = 0;
	adjMX[10][10] = 0;
	adjMX[10][11] = 11;

	adjMX[11][0] = 0;
	adjMX[11][1] = 0;
	adjMX[11][2] = 0;
	adjMX[11][3] = 0;
	adjMX[11][4] = 0;
	adjMX[11][5] = 0;
	adjMX[11][6] = 0;
	adjMX[11][7] = 0;
	adjMX[11][8] = 0;
	adjMX[11][9] = 0;
	adjMX[11][10] = 0;
	adjMX[11][11] = 0;

	initial = 0;
	end = 11;
	Dijkstra(adjMX, n, initial, end);

	return 0;
	}
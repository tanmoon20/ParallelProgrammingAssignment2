#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <map>

using namespace std;

struct coordinate {
	string name = "";
	int x = 0;
	int y = 0;
};

struct vertexEnd
{
	string destinationName;
	double distance = -1.0;
	map<int, bool> extinguished;
	vector<coordinate> fireExtinguished;
};

struct vertexStart {
	string startName;
	vector<vertexEnd> edgeList;
};

struct route {
	vector<int> passedVertex;
	vector<int> extinguishedFire;
	map<int, bool> extinguished;
};

struct FireExtinguished
{
	map<int, bool> extinguished;
	int totalExtinguished = 0;
};

coordinate airport;
vector<coordinate> fire;
vector<coordinate> shortest;
double minDistance = INFINITY;
vector<vertexStart> distanceGraph;
//-------------new!!!
//store coordinate and name of terrain points
vector<coordinate> terrain;
//store the number of fire name;
int* shortestAirplaneToAirplane;
int* fires;
route** shortestRoutes;
double** distanceMatrix;
map<int, bool>** extinguishedMatrix;
//record fireExtinguished of current permutation of fires
FireExtinguished extinguished;

double getTotalDistance(vector<coordinate>&);
void extinguishFire(int, int, vector<coordinate>&, int);
double getDistanceFromDistanceTable(int, int);
void heapPerm(int);
void readLocation(void);
void saveLocation(void);
double getDistanceBetweenPoints(coordinate, coordinate);
double getDistanceBetweenPointsNew(int start, int end);
double getShortestDistanceFromPointToLine(double, double, double, coordinate);
int getVertexIndex(string);
void exchangeValue(int&, int&);
bool definitelyLessThan(float, float);
bool essentiallyEqual(float, float);
double getAbsolute(double x);
void createExtinguishTable(int, int);
void createExtinguishTableNew();
void heapPermNew(int length);
double getTotalDistanceNew(void);
bool extinguishFireNew(int indexStart, int indexEnd, int firePoint);

bool extinguishFireNew(int indexStart, int indexEnd, int firePoint) {

	return extinguishedMatrix[indexStart][indexEnd].find(firePoint)->second;
}

double getTotalDistanceNew(void) {
	int firesSize = terrain.size() - 1;
	double totalDistance = 0;
	//get the distance between airport and the first fire
	totalDistance += distanceMatrix[0][fires[0]];

	for (int i = 1; i <= firesSize; i++)
	{
		extinguished.extinguished.insert({ i,false });
	}
	

	//if terrain has more than 1 fire
	if (terrain.size() > 2) {
		//extinguishFire between airport and first fire
		for (int i = 0; i < firesSize; i++)
		{
			if (extinguishedMatrix[0][fires[0]].find(fires[i])->second)
			{
				extinguished.extinguished.find(fires[i])->second = true;
				extinguished.totalExtinguished++;
			}
		}
	}
	else
		return totalDistance;

	vector<coordinate> fireTemp2;
	int indexPreviousFireTemp;
	int indexCurrentFireTemp;

	// extinguish fire between fire points
	for (int i = 1; i < firesSize; i++)
	{
		//visit fire not extinguished
		if (!extinguished.extinguished.find(fires[i])->second) {
			totalDistance += distanceMatrix[fires[i - 1]][fires[i]];
			//check if fire ahead is extinguishable
			if (i + 1 < firesSize)
			{
				//stop calculate if larger than minDist	
				if (totalDistance >= minDistance)
				{
					return totalDistance;
				}

				if (extinguishFireNew(i - 1, i, i + 1))
				{
					extinguished.extinguished.find(fires[i])->second = true;
					extinguished.totalExtinguished++;
				}
				
				//check if the point can reach to the end(airport) and extinguish all remaining fire
				if (extinguished.totalExtinguished > 2)
				{
					for (int j = i + 1; j < firesSize; j++)
					{
						if (!extinguished.extinguished.find(fires[i])->second)
						{
							if (extinguishFireNew(0, i, j))
							{
								if (j == firesSize - 1)
								{
									return totalDistance + distanceMatrix[0][fires[i]];
								}
							}
							else
								break;
						}
					}
				}
			}
		}
	}

	//get the distance between airport and the last fire
	totalDistance += distanceMatrix[0][fires[firesSize-1]];
	return totalDistance;
}


void heapPermNew(int length)
{
	if (length == 1)
	{
		for (int i = 0; i < terrain.size() - 2; i++)
		{
			cout << fires[i] << " ";
		}
		cout << endl;
		//compute shortest path that can extinguish all fire
		double currentDistance = getTotalDistanceNew();
		if (currentDistance < minDistance) {
			minDistance = currentDistance;
			for (int i = 1; i < terrain.size(); i++)
			{
				if (!extinguished.extinguished.find(i - 1)->second)
				{
					shortestAirplaneToAirplane[i] = 1;
				}
				else
				{
					shortestAirplaneToAirplane[i] = 0;
				}
			}
		}
		extinguished.extinguished.clear();
		extinguished.totalExtinguished = 0;
	}
	else
	{
		length -= 1;
		heapPermNew(length);
		for (int i = 0; i < length; i++) {
			if (length % 2 != 0)
			{
				swap(fires[i], fires[length]);
			}
			else
			{
				swap(fires[0], fires[length]);
			}
			heapPermNew(length);
		}
	}
}

void createExtinguishTableNew(void) {
	int size = terrain.size();
	extinguishedMatrix = new map<int, bool>* [size];

	map<int, bool> extinguished;
	for (int i = 0; i < size; i++) {
		extinguished.insert({ i,false });
	}

	for (int i = 0; i < size; i++) {
		extinguishedMatrix[i] = new map<int, bool> [size];
		for (int j = 0; j < size; j++) {
			extinguishedMatrix[i][j] = extinguished;
		}
	}

	for (int i = 1; i < terrain.size(); i++) 
	{
		coordinate start = terrain.at(i);
		coordinate end = terrain.at(i - 1);

		int pointStart = i - 1;
		int pointEnd = i;

		//get the gradient of start -> end line
		double gradientLine = (start.y - end.y) * 1.0 / (start.x - end.x);

		//find shortest distance between fire point and start->end line
		//get constant to get eqn of start -> end line (y = gradientLine(x) + constantLine)
		double constantLine = start.y - (gradientLine * start.x);

		coordinate tempFirePointCoordinate;
		//loop all the fire around to check if extinguishable
		for (int j = 1; j < terrain.size(); j++)
		{
			if (i != pointStart && i != pointEnd)
			{
				tempFirePointCoordinate = terrain.at(i);
				//calculate shortest distance between fire point to start->end line
				double shortestDistFireToLine = getShortestDistanceFromPointToLine(gradientLine, -1, constantLine, tempFirePointCoordinate);
				//check if the y-distance between point and start->end line < 50, if yes, extinguish fire
				if (shortestDistFireToLine <= 50)
				{
					//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
					//get eqn of left, right line that are perpendicular to start->end line and intersects start/end point.
					double gradientLeftRight = (-1.0 / gradientLine);
					double constantStart = start.y - (gradientLeftRight * start.x);
					double constantEnd = end.y - (gradientLeftRight * end.x);

					//calculate shortest distance between fire point to left/right line
					double shortestDistPointToStartLine = getShortestDistanceFromPointToLine(gradientLeftRight, -1, constantStart, tempFirePointCoordinate);
					double shortestDistPointToEndLine = getShortestDistanceFromPointToLine(gradientLeftRight, -1, constantEnd, tempFirePointCoordinate);
					double sumOfshortestDistPointToStartAndEndLine = shortestDistPointToStartLine + shortestDistPointToEndLine;
					//the index passed is the index terrain of distanceGraph
					double distanceBetweenPoints = distanceMatrix[pointStart][pointEnd];

					//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
					if (definitelyLessThan(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints) || essentiallyEqual(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints))
					{
						//record fire extinguished into table
						extinguishedMatrix[pointStart][pointEnd].find(i)->second = true;
						//distanceGraph.at(pointStart).edgeList.at(pointEnd).extinguished.find(i)->second = true;
					}
				}
			}
		}
	}

}

double getDistanceBetweenPointsNew(int start, int end)
{
	double x = terrain.at(end).x - terrain.at(start).x;
	double y = terrain.at(end).y - terrain.at(start).y;
	return sqrt((x * x) + (y * y));
}

void performFloydWarshall() {
	int size = terrain.size();
	distanceMatrix = new double* [size];
	shortestRoutes = new route * [size];


	for (int i = 0; i < size; i++) {
		distanceMatrix[i] = new double[size];
		shortestRoutes[i] = new route[size];
	}

	map<int, bool> extinguished;
	for (int i = 0; i < size; i++) {
		extinguished.insert({ i,false });
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			shortestRoutes[i][j].extinguished = extinguished;
			if (i == j) {
				distanceMatrix[i][j] = 0;
			}
			else
			{
				distanceMatrix[i][j] = getDistanceBetweenPointsNew(i,j);
			}
		}
	}

	//int cardinality = size;
	//vector<coordinate> temp;
	//temp.reserve(3);
	//for (int k = 0; k < cardinality; k++) 
	//{
	//	for (int i = 0; i < cardinality; i++)
	//	{
	//		for (int j = 0; j < cardinality; j++)
	//		{
	//			cout << i << " " << j << " " << k << endl;
	//			if (i != j) {
	//				
	//				bool swapped = false;
	//				if (i > j) {
	//					swap(i, j);
	//					swapped = true;
	//				}
	//				cout << i << " == " << j << " == " << k << endl;
	//				if (distanceGraph.at(i).edgeList.at(j - (i + 1)).extinguished.find(k)->second)
	//				{
	//					cout << "yes" << endl;
	//					shortestRoutes[i][j].extinguished.find(k)->second = true;
	//					shortestRoutes[i][j].extinguishedFire.push_back(k);
	//				}
	//				if(swapped)
	//					swap(i, j);
	//				cout << i << " == " << j << " == " << k << endl;
	//			}

	//			/*for (int l = 0; l < distanceGraph.at(i).edgeList.at(j - (i + 1)).fireExtinguished.size(); l++) {
	//				if (distanceGraph.at(i).edgeList.at(j - (i + 1)).fireExtinguished.at(l).name.compare(fire.at(k - 1).name) == 0)
	//				{
	//					shortestRoutes[i][j].extinguishedFire.push_back(k);

	//				}
	//				else
	//				{
	//					if (distance[i][j] > distance[i][k] + distance[k][j])
	//					{
	//						shortestRoutes[i][j].passedVertex.push_back(k);
	//					}
	//				}
	//				distance[i][j] = distance[i][k] + distance[k][j];
	//			}*/
	//			/*if (distance[i][j] > distance[i][k] + distance[k][j]) {
	//				distance[i][j] = distance[i][k] + distance[k][j];
	//				cout << "aaa";
	//			}*/

	//		}
	//	}
	//}

	//free memory
	//for (int i = 0; i < size; i++)
	//{
	//	delete[] distanceMatrix[i];
	//}
	//delete[] distanceMatrix;
	delete[] shortestRoutes;
}

int main(void) {

	readLocation();
	performFloydWarshall();
	createExtinguishTableNew();
	//there exists fire
	if (terrain.size() > 1) {
		shortestAirplaneToAirplane = new int(terrain.size());
		shortestAirplaneToAirplane[0] = 1;
		heapPermNew(terrain.size() - 1);
	}
	saveLocation();
	
	return 0;
}

double getTotalDistance(vector<coordinate>& fireTemp) {
	double totalDistance = 0;
	int indexFirstFireTemp = getVertexIndex(fireTemp.front().name);
	int indexLastFireTemp;
	int indexAirport = 0;
	//get the distance between airport and the first fire
	totalDistance += getDistanceFromDistanceTable(indexAirport, indexFirstFireTemp);

	if (fireTemp.size() > 1) {
		//extinguishFire between airport and first fire
		extinguishFire(indexAirport, indexFirstFireTemp, fireTemp, 1);
	}
	//only a fire in the terrain
	else
		return totalDistance;

	vector<coordinate> fireTemp2;
	int indexPreviousFireTemp;
	int indexCurrentFireTemp;

	// extinguish fire between fire points
	for (int i = 1; i < fireTemp.size(); i++)
	{
		indexPreviousFireTemp = getVertexIndex(fireTemp.at(i - 1).name);
		indexCurrentFireTemp = getVertexIndex(fireTemp.at(i).name);
		totalDistance += getDistanceFromDistanceTable(indexPreviousFireTemp, indexCurrentFireTemp);
		if (i + 1 < fireTemp.size())
		{
			extinguishFire(indexPreviousFireTemp, indexCurrentFireTemp, fireTemp, i + 1);
			//stop calculate if larger than minDist	
			if (totalDistance >= minDistance)
			{
				return totalDistance;
			}
			//check if the point can reach to the end(airport) and extinguish all remaining fire
			if (fireTemp.size() > 2)
			{
				//the first point ady calculated with end(start(airport)), so calculate starting from 2nd point, i = 1
				fireTemp2 = fireTemp;
				extinguishFire(indexAirport, getVertexIndex(fireTemp2.at(i).name), fireTemp2, i + 1);
				//if all fire can be extinguished, go back to airport
				if (fireTemp2.size() == i + 1)
				{
					fireTemp = fireTemp2;
					indexLastFireTemp = getVertexIndex(fireTemp.back().name);
					totalDistance += getDistanceFromDistanceTable(indexAirport, indexLastFireTemp);
					return totalDistance;
				}
				else
					fireTemp2 = fireTemp;
			}
		}
	}

	indexLastFireTemp = getVertexIndex(fireTemp.back().name);
	//get the distance between airport and the last fire
	totalDistance += getDistanceFromDistanceTable(indexAirport, indexLastFireTemp);
	return totalDistance;
}

void extinguishFire(int indexStart, int indexEnd, vector<coordinate>& fireTemp, int surroundingFireIndex) {

	exchangeValue(indexStart, indexEnd);
	indexEnd -= (indexStart + 1);

	for (int i = 0; i < distanceGraph.at(indexStart).edgeList.at(indexEnd).fireExtinguished.size(); i++)
	{
		for (int j = surroundingFireIndex; j < fireTemp.size(); j++)
		{
			if (distanceGraph.at(indexStart).edgeList.at(indexEnd).fireExtinguished.at(i).name.compare(fireTemp.at(j).name) == 0)
			{
				//fire extinguished is found
				surroundingFireIndex = j;
				fireTemp.erase(fireTemp.begin() + j);
				//escape for loop
				break;
			}
		}
	}
}

double getDistanceFromDistanceTable(int indexStart, int indexEnd) {
	exchangeValue(indexStart, indexEnd);
	//change the index of distanceGraph to the index of edgeList
	indexEnd -= (indexStart + 1);
	double distanceBetweenPoints = distanceGraph.at(indexStart).edgeList.at(indexEnd).distance;
	return distanceBetweenPoints;
}

void heapPerm(int length)
{
	if (length == 1)
	{
		vector<coordinate> temp;
		temp = fire;
		//compute shortest path that can extinguish all fire
		double currentDistance = getTotalDistance(temp);
		if (currentDistance < minDistance) {
			minDistance = currentDistance;
			shortest = temp;
		}

	}
	else
	{
		coordinate dummy;
		length -= 1;
		heapPerm(length);
		for (int i = 0; i < length; i++) {
			dummy = fire.at(length);
			if (length % 2 != 0)
			{
				swap(fire.at(i), fire.at(length));
			}
			else
			{
				swap(fire.front(), fire.at(length));
			}
			heapPerm(length);
		}
	}
}

void readLocation(void) {

	string txt;
	ifstream inFile("terrain.txt");
	coordinate terrainCoordinate;
	if (!inFile.is_open()) {
		cout << "Terrain file is not found!.\n";
		return;
	}
	else
	{
		while (inFile >> txt)
		{
			//is tree
			if (txt[0] == 'T')
			{
				inFile >> txt;
				inFile >> txt;
				continue;
			}
			//is airport/fire
			else
			{
				terrainCoordinate.name = txt;
				inFile >> terrainCoordinate.x;
				inFile >> terrainCoordinate.y;
			}
			//create terrain
			terrain.push_back(terrainCoordinate);
		}
		inFile.close();
	}

	//create fires
	fires = new int(terrain.size() - 1);
	for (int i = 1; i < terrain.size(); i++)
	{
		fires[i - 1] = i;
	}

}

void saveLocation(void) {
	ofstream outFile("solution.txt");
	outFile << airport.name << endl;
	for (unsigned int i = 0; i < shortest.size(); i++)
	{
		outFile << shortest.at(i).name << endl;
	}
	outFile << airport.name << endl;
	outFile.close();
}

double getDistanceBetweenPoints(coordinate start, coordinate end)
{
	double x = end.x - start.x;
	double y = end.y - start.y;
	return sqrt((x * x) + (y * y));
}

//use formula of shortest distance from point to line
double getShortestDistanceFromPointToLine(double A, double B, double C, coordinate point)
{
	double shortestDistPartA = (A * point.x) + (B * point.y) + C;
	double shortestDistPartB = getAbsolute(shortestDistPartA) / (sqrt((A * A) + (B * B)));
	return shortestDistPartB;
}

//get index of vertex in the vector
int getVertexIndex(string name)
{
	//is airport, ascii of t is 116
	if (!name.compare(airport.name))
	{
		return 0;
	}
	else
	{
		for (int i = 1; i < distanceGraph.size(); i++)
		{
			if (!distanceGraph.at(i).startName.compare(name))
			{
				return i;
			}
		}
	}
}

//index with smaller value will be startVertex
void exchangeValue(int& startVertex, int& endVertex)
{
	int dummy;
	if (startVertex > endVertex)
	{
		dummy = startVertex + endVertex;
		startVertex = endVertex;
		endVertex = dummy - startVertex;
	}
}

//compare double to check less than
bool definitelyLessThan(float a, float b)
{
	float absA = getAbsolute(a);
	float absB = getAbsolute(b);
	bool isLessThan = (b - a) > ((absA < absB ? absB : absA) * numeric_limits<double>::epsilon());
	return isLessThan;
}

//compare double to check equal
bool essentiallyEqual(float a, float b)
{
	float absA = getAbsolute(a);
	float absB = getAbsolute(b);
	bool isEqual = getAbsolute(a - b) <= ((absA > absB ? absB : absA) * numeric_limits<double>::epsilon());
	return isEqual;
}

double getAbsolute(double x)
{
	if (x < -1)
		x *= -1.0;
	return x;
}

void createExtinguishTable(int indexStart, int indexEnd) {
	coordinate start, end;
	if (indexStart == 0)
	{
		start = airport;
	}
	else {
		start = fire.at(indexStart - 1);
	}
	//current indexEnd value is the index in distanceGraph
	end = fire.at(indexEnd - 1);
	//change indexEnd to index of edgeList
	indexEnd -= (indexStart + 1);

	//get the gradient of start -> end line
	double gradientLine = (start.y - end.y) * 1.0 / (start.x - end.x);

	//find shortest distance between fire point and start->end line
	//get constant to get eqn of start -> end line (y = gradientLine(x) + constantLine)
	double constantLine = start.y - (gradientLine * start.x);

	coordinate tempFirePointCoordinate;
	for (int i = 0; i < fire.size(); i++)
	{
		if (!(fire.at(i).name.compare(start.name) == 0 || fire.at(i).name.compare(end.name) == 0))
		{
			tempFirePointCoordinate = fire.at(i);
			//calculate shortest distance between fire point to start->end line
			double shortestDistFireToLine = getShortestDistanceFromPointToLine(gradientLine, -1, constantLine, tempFirePointCoordinate);
			//check if the y-distance between point and start->end line < 50, if yes, extinguish fire
			if (shortestDistFireToLine <= 50)
			{
				//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
				//get eqn of left, right line that are perpendicular to start->end line and intersects start/end point.
				double gradientLeftRight = (-1.0 / gradientLine);
				double constantStart = start.y - (gradientLeftRight * start.x);
				double constantEnd = end.y - (gradientLeftRight * end.x);

				//calculate shortest distance between fire point to left/right line
				double shortestDistPointToStartLine = getShortestDistanceFromPointToLine(gradientLeftRight, -1, constantStart, tempFirePointCoordinate);
				double shortestDistPointToEndLine = getShortestDistanceFromPointToLine(gradientLeftRight, -1, constantEnd, tempFirePointCoordinate);
				double sumOfshortestDistPointToStartAndEndLine = shortestDistPointToStartLine + shortestDistPointToEndLine;
				//the index passed is the index terrain of distanceGraph
				double distanceBetweenPoints = getDistanceFromDistanceTable(indexStart, indexEnd + indexStart + 1);

				//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
				if (definitelyLessThan(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints) || essentiallyEqual(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints))
				{
					//record fire extinguished into table
					//the indexStart is index of distanceGraph, indexEnd is index of edgeList;
					distanceGraph.at(indexStart).edgeList.at(indexEnd).fireExtinguished.push_back(tempFirePointCoordinate);
					distanceGraph.at(indexStart).edgeList.at(indexEnd).extinguished.find(i + 1)->second = true;
				}
			}
		}
	}
}
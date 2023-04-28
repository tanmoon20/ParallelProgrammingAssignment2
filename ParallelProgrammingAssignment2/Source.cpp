#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>

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
	vector<coordinate> fireExtinguished;
};

struct vertexStart {
	string startName;
	vector<vertexEnd> edgeList;
};

coordinate airport;
vector<coordinate> fire;
vector<coordinate> shortest;
double minDistance = INT64_MAX;
vector<vertexStart> distanceGraph;

double getTotalDistance(vector<coordinate>&);
void extinguishFire(int, int, vector<coordinate>&, int);
double getDistanceFromDistanceTable(int, int);
void heapPerm(int, vector<coordinate>&);
bool readLocation(void);
void saveLocation(void);
double getDistanceBetweenPoints(coordinate, coordinate);
double getShortestDistanceFromPointToLine(double, double, double, coordinate);
int getVertexIndex(string);
void exchangeValue(int&, int&);
bool definitelyLessThan(float, float);
bool essentiallyEqual(float, float);
double getAbsolute(double x);
void createExtinguishTable();
void createDistanceTable();

int main(void) {
	auto start = chrono::high_resolution_clock::now();

	if (readLocation())
	{
		createDistanceTable();
		createExtinguishTable();

		//there exists fire
		if (distanceGraph.size() > 1) {
			//parallel permutation to do brute force
#pragma omp parallel for
			for (int i = 0; i < fire.size(); i++) {
				vector<coordinate> firedupe = fire;
				swap(firedupe.at(i), firedupe.back());
				heapPerm(firedupe.size() - 1, firedupe);
			}
		}

		saveLocation();
	}

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
		//stop calculate if larger than minDist	
		if (totalDistance >= minDistance)
		{
			return totalDistance;
		}

		if (i + 1 < fireTemp.size())
		{
			extinguishFire(indexPreviousFireTemp, indexCurrentFireTemp, fireTemp, i + 1);
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

void heapPerm(int length, vector<coordinate>& toPermute)
{
	if (length == 1)
	{
		vector<coordinate> temp;
		temp = toPermute;

		//compute shortest path that can extinguish all fire
		double currentDistance = getTotalDistance(temp);
#pragma omp critical
		if (currentDistance < minDistance) {
			minDistance = currentDistance;
			shortest = temp;
		}

	}
	else
	{
		length -= 1;
		heapPerm(length, toPermute);
		for (int i = 0; i < length; i++) {
			if (length % 2 != 0)
			{
				swap(toPermute.at(i), toPermute.at(length));
			}
			else
			{
				swap(toPermute.front(), toPermute.at(length));
			}
			heapPerm(length, toPermute);
		}
	}
}

bool readLocation(void) {
	struct coordinate temp;
	vertexStart pointInTerrain;
	string txt;
	ifstream inFile("terrain.txt");

	if (!inFile.is_open()) {
		cout << "Terrain file is not found!.\n";
		return 0;
	}
	else
	{
		while (inFile >> txt)
		{
			//is fire
			if (txt[0] == 'F')
			{
				temp.name = txt;
				inFile >> temp.x;
				inFile >> temp.y;
				fire.push_back(temp);

				pointInTerrain.startName = temp.name;
			}
			//is tree
			else if (txt[0] == 'T')
			{
				inFile >> txt;
				inFile >> txt;
				continue;
			}
			//is airport
			else
			{
				airport.name = txt;
				inFile >> airport.x;
				inFile >> airport.y;

				pointInTerrain.startName = airport.name;
			}
			distanceGraph.push_back(pointInTerrain);
		}
		inFile.close();
	}
	return 1;
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

void createExtinguishTable() {

	for (int iStart = 0; iStart < distanceGraph.size(); iStart++)
	{
		for (int iEnd = iStart + 1; iEnd < distanceGraph.size(); iEnd++)
		{
			int indexStart = iStart;
			int indexEnd = iEnd;
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
						}
					}
				}
			}
		}
	}
}

void createDistanceTable() {
	//insert all vertex name into distanceGraph
	vertexEnd destination;

	for (int i = 0; i < distanceGraph.size(); i++)
	{
		distanceGraph.at(i).edgeList.reserve(distanceGraph.size() - (i + 1));
		for (int j = i + 1; j < distanceGraph.size(); j++)
		{
			//create distance table between two points of all points
			if (i == 0) {
				destination.distance = getDistanceBetweenPoints(airport, fire.at(j - 1));
			}
			else {
				destination.distance = getDistanceBetweenPoints(fire.at(i - 1), fire.at(j - 1));
			}
			destination.destinationName = distanceGraph.at(j).startName;
			distanceGraph.at(i).edgeList.push_back(destination);
		}

	}
}
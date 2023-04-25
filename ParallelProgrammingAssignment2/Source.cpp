#include <iostream>
#include <list>
#include <fstream>
#include <vector>

using namespace std;

struct coordinate {
	string name = "";
	int x = 0;
	int y = 0;
}airport;

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

vector<coordinate> terrainCoordinate;
vector<coordinate> fire;
vector<coordinate> shortest;
double minDistance = INT64_MAX;
vector<vertexStart> distancetable;

struct edge;
struct vertex;

struct edge {
	int val;
	double cost = 0;
	list<int> fireExtinguished;
	list<int> shortestPath;
};

//store adjacency list
struct vertex {
	//val = index, 0 = airport, 1 = Fire1, 2 = Fire2 etc
	int val;
	list<edge> edgeList;
};

class Graph {

public:
	list<vertex> vertexList;
	int size() {
		return vertexList.size();
	}

	void insertVertex(int val) {
		vertex v;
		v.val = val;
		vertexList.push_back(v);
	}

	void insertEdge(int vertexVal, int edgeVal) {
		vertex v;
		edge e;

		e.val = edgeVal;
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				i.edgeList.push_back(e);
				return;
			}
		}
	}

	void insertEdge(int vertexVal, int edgeVal, double distance) {
		edge e;

		e.val = edgeVal;
		e.cost = distance;
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				i.edgeList.push_back(e);
				return;
			}
		}
	}

	double getDistance(int vertexVal, int edgeVal) {
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					if (j.val == edgeVal)
						return j.cost;
				}
			}
		}
	}

	void setDistance(int vertexVal, int edgeVal, double cost) {
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					if (j.val == edgeVal)
						j.cost = cost;
				}
			}
		}
	}

	list<int> getShortestPath(int vertexVal, int edgeVal) {
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					if (j.val == edgeVal)
						return j.shortestPath;
				}
			}
		}
	}

	void setFireExtinguished(int vertexVal, int edgeVal, int fireExtinguished)
	{
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					/*if (fireExtinguished > j.fireExtinguished.back())
					{*/
					j.fireExtinguished.push_back(fireExtinguished);
					//j.fireExtinguished.sort();
					return;
					//}
				}
			}
		}
	}

	void setFireExtinguished(int vertexVal, int edgeVal, list<int> fireExtinguished)
	{
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					j.fireExtinguished = fireExtinguished;
					//j.fireExtinguished.sort();
					return;
				}
			}
		}
	}

	int getFireExtinguishedMaxValue(int vertexVal, int edgeVal) {
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					return j.fireExtinguished.back();
				}
			}
		}
	}

	list<int> getFireExtinguished(int vertexVal, int edgeVal) {
		for (auto& i : vertexList) {
			if (i.val == vertexVal) {
				for (auto& j : i.edgeList)
				{
					return j.fireExtinguished;
				}
			}
		}
	}
};

Graph terrainMap;
Graph distanceGraph;

double getTotalDistance(vector<coordinate>&);
void extinguishFire(int, int, vector<coordinate>&, int);
double getDistanceFromDistanceTable(int, int);
void heapPerm(int);
void readLocation(void);
void saveLocation(void);
double getDistanceBetweenPoints(coordinate, coordinate);
double getShortestDistanceFromPointToLine(double, double, double, coordinate);
int getVertexIndex(string);
void exchangeValue(int&, int&);
bool definitelyLessThan(float, float);
bool essentiallyEqual(float, float);
double getAbsolute(double x);
void createExtinguishTable(int, int);
void floydWarshallAlgo();
void createExtinguishTablev2(int indexStart, int indexEnd);
double getDistanceFromDistanceTablev2(int indexStart, int indexEnd);

void readLocation(void) {
	struct coordinate temp;
	vertexStart pointInTerrain;

	string txt;
	ifstream inFile("terrain.txt");
	if (!inFile.is_open()) {
		cout << "Terrain file is not found!.\n";
		return;
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
				terrainCoordinate.push_back(temp);
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
				terrainCoordinate.push_back(airport);
			}
			distancetable.push_back(pointInTerrain);
		}
		inFile.close();
	}

	//insert all vertex name into distanceGraph (X check)
	vertexEnd destination;
	for (int i = 0; i < distancetable.size(); i++)
	{
		distancetable.at(i).edgeList.reserve(distancetable.size() - (i + 1));
		for (int j = i + 1; j < distancetable.size(); j++)
		{
			//create distance table between two points of all points
			if (i == 0) {
				destination.distance = getDistanceBetweenPoints(airport, fire.at(j - 1));
			}
			else {
				destination.distance = getDistanceBetweenPoints(fire.at(i - 1), fire.at(j - 1));
			}
			destination.destinationName = distancetable.at(j).startName;
			distancetable.at(i).edgeList.push_back(destination);
		}

	}

	//insert vertex, edge and weight to distanceGraph
	for (int i = 0; i < terrainCoordinate.size(); i++)
	{
		distanceGraph.insertVertex(i);
		for (int j = 0; j < terrainCoordinate.size(); j++)
		{
			double distance;
			//when vertex = edge, the distance = 0 by default
			if (i == j) {
				distanceGraph.insertEdge(i, j);
			}
			else if (i < j) {
				distance = getDistanceBetweenPoints(terrainCoordinate.at(i), terrainCoordinate.at(j));
				//createExtinguishTablev2(i, j);
				distanceGraph.insertEdge(i, j, distance);
			}
			else
			{
				distance = distanceGraph.getDistance(j, i);
				//distanceGraph.setFireExtinguished(j, i, distanceGraph.getFireExtinguished(j,i));
				distanceGraph.insertEdge(i, j, distance);
			}
		}
	}

	//insert fireExtinguished to distanceGraph
	for (int i = 0; i < terrainCoordinate.size(); i++)
	{
		for (int j = 0; j < terrainCoordinate.size(); j++)
		{
			if (i < j) {
				createExtinguishTablev2(i, j);
			}
			else if (i > j)
			{
				distanceGraph.setFireExtinguished(j, i, distanceGraph.getFireExtinguished(j, i));
			}
		}
	}


	//floydWarshallAlgo();

	//create extinguish table
	for (int i = 1; i < distancetable.size(); i++)
	{
		for (int j = i + 1; j < distancetable.size(); j++)
		{
			//create distance table between two points of all points
			createExtinguishTable(i - 1, j);
		}

	}
}

void createExtinguishTablev2(int indexStart, int indexEnd) {
	coordinate start, end;
	start = terrainCoordinate.at(indexStart);
	end = terrainCoordinate.at(indexEnd);

	//get the gradient of start -> end line
	double gradientLine = (start.y - end.y) * 1.0 / (start.x - end.x);

	//find shortest distance between fire point and start->end line
	//get constant to get eqn of start -> end line (y = gradientLine(x) + constantLine)
	double constantLine = start.y - (gradientLine * start.x);

	coordinate tempFirePointCoordinate;
	//starting from firePoint
	for (int i = 1; i < terrainCoordinate.size(); i++)
	{
		if (!(i == indexStart || i == indexEnd))
		{
			tempFirePointCoordinate = terrainCoordinate.at(i);
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
				double distanceBetweenPoints = getDistanceFromDistanceTablev2(indexStart, indexEnd);
				cout << distanceBetweenPoints << endl;
				//double distanceBetweenPoints = getDistanceFromDistanceTable(indexStart, indexStart + 1);

				//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
				if (definitelyLessThan(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints) || essentiallyEqual(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints))
				{
					//record fire extinguished into table
					//the indexStart is index of distanceGraph, indexEnd is index of edgeList;
					//distancetable.at(indexStart).edgeList.at(indexEnd).fireExtinguished.push_back(tempFirePointCoordinate);
					//firei = i
					distanceGraph.setFireExtinguished(indexStart, indexEnd, i);
				}
			}
		}
	}
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
	//indexEnd is the index in edgeList
	end = fire.at(indexEnd - 1);
	//change indexEnd to index of EdgeList
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
				//double distanceBetweenPoints = getDistanceFromDistanceTable(indexStart, indexStart + 1);

				//check if the x-coordinate of point is within start->end line, if yes, extinguish fire
				if (definitelyLessThan(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints) || essentiallyEqual(sumOfshortestDistPointToStartAndEndLine, distanceBetweenPoints))
				{
					//record fire extinguished into table
					//the indexStart is index of distanceGraph, indexEnd is index of edgeList;
					distancetable.at(indexStart).edgeList.at(indexEnd).fireExtinguished.push_back(tempFirePointCoordinate);
					//fire1 = i + 1(0+1)
					//distanceGraph.setFireExtinguished(indexStart, indexEnd + (indexStart + 1), i + 1);
				}
			}
		}
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
	return x;
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
		for (int i = 1; i < distancetable.size(); i++)
		{
			if (!distancetable.at(i).startName.compare(name))
			{
				return i;
			}
		}
	}
}

//index with smaller value will be startVertex
void exchangeValue(int& startVertex, int& endVertex)
{
	if (startVertex > endVertex)
	{
		swap(startVertex, endVertex);
	}
}

double getDistanceFromDistanceTable(int indexStart, int indexEnd) {
	exchangeValue(indexStart, indexEnd);
	//change the index of distanceGraph to the index of edgeList
	indexEnd -= (indexStart + 1);
	double distanceBetweenPoints = distancetable.at(indexStart).edgeList.at(indexEnd).distance;
	return distanceBetweenPoints;
}

double getDistanceFromDistanceTablev2(int indexStart, int indexEnd) {
	//change the index of distanceGraph to the index of edgeList
	double distanceBetweenPoints = distanceGraph.getDistance(indexStart, indexEnd);
	cout << distanceBetweenPoints << endl;
	return distanceBetweenPoints;
}

void floydWarshallAlgo() {
	cout << "a";
	for (int i = 0; i < distanceGraph.size(); i++)
	{
		for (int j = 0; j < distanceGraph.size(); j++)
		{
			cout << "b";
			double distance;
			if (i == j) {
				distance = 0;
			}
			else if (i > j) {
				distance = distancetable.at(j).edgeList.at(i - j - 1).distance;
			}
			else
			{
				distance = distancetable.at(i).edgeList.at(j - i - 1).distance;
			}
			distanceGraph.insertEdge(i, j, distance);
		}
	}

	list<int> shortestPath;
	list<int> fireExtinguished;

	for (int k = 0; k < distanceGraph.size(); k++)
	{
		for (int i = 0; i < distanceGraph.size(); i++)
		{
			for (int j = 0; j < distanceGraph.size(); j++)
			{
				if (distanceGraph.getDistance(i, j) > distanceGraph.getDistance(i, k) + distanceGraph.getDistance(k, j))
				{

					distanceGraph.setDistance(i, j, distanceGraph.getDistance(i, k) + distanceGraph.getDistance(k, j));
					//shortestPath.push_back(k);
				}
			}
		}
	}
}

int main() {
	readLocation();

	//for (int i = 0; i < distancetable.size(); i++)
	//{
	//	terrainMap.insertVertex(i);
	//	distanceGraph.insertVertex(i);
	//	for (int j = 0; j < distancetable.at(i).edgeList.size();j++)
	//	{
	//		//airport -> edgeList store index 1 to ...
	//		//fire1 -> edgeList store index 2 to ...
	//		terrainMap.insertEdge(i, i + j + 1, distancetable.at(i).edgeList.at(j).distance);
	//	}
	//}




	return 0;
}
//------------------------------------------------------------------------------
// exampleapp.cc
// (C) 2015-2017 Individual contributors, see AUTHORS file
//------------------------------------------------------------------------------
#include "config.h"
#include "exampleapp.h"
#include <cstring>
#include <glm/glm.hpp>
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <ctime>
#include "imgui.h"

#define STRING_BUFFER_SIZE 8192
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <vector>
const GLchar* vs =
"#version 430\n"
"layout(location=0) in vec3 pos;\n"
"layout(location=1) in vec4 color;\n"
"layout(location=0) out vec4 Color;\n"
//"uniform mat4 rotationMtrx;\n"
"void main()\n"
"{\n"
//"   gl_Position = rotationMtrx * vec4(pos, 1.0);\n"
"	gl_Position = vec4(pos.x, pos.y, pos.z, 1);\n"
"	Color = color;\n"
"}\n";

const GLchar* ps =
"#version 430\n"
"layout(location=0) in vec4 color;\n"
"out vec4 Color;\n"
"void main()\n"
"{\n"
"	Color = color;\n"
"}\n";


using namespace Display;
namespace Example
{

	//------------------------------------------------------------------------------
	/**
	*/
	ExampleApp::ExampleApp()
	{
		// empty
	}

	//------------------------------------------------------------------------------
	/**
	*/
	ExampleApp::~ExampleApp()
	{
		// empty
	}

	//------------------------------------------------------------------------------
	/**
	*/


	int rec = 1;
	std::vector<GLfloat> baseTri;


	std::vector<GLfloat> getMahPointsoutline(glm::vec2 startPoint, glm::vec2 endPoint, glm::vec2 basePoint, int rec) {
		glm::vec2 firstThird = { ((2 * startPoint.x) + endPoint.x) / 3, ((2 * startPoint.y) + endPoint.y) / 3 };
		glm::vec2 secondThird = { ((2 * endPoint.x) + startPoint.x) / 3, ((2 * endPoint.y) + startPoint.y) / 3 };
		GLfloat Ux = endPoint.x - startPoint.x;
		GLfloat Uy = endPoint.y - startPoint.y;
		GLfloat Vx = endPoint.x - startPoint.x;
		GLfloat Vy = startPoint.y - endPoint.y;
		GLfloat Qx = startPoint.x + (1.0f / 2.0f)*Ux + (sqrt(3.0f) / 6.0f)*Vy;
		GLfloat Qy = startPoint.y + (1.0f / 2.0f)*Uy + (sqrt(3.0f) / 5.0f)*Vx;
		glm::vec2 midPoint = { Qx, Qy };
		glm::vec2 base1 = { ((2 * endPoint.x) + basePoint.x) / 3, ((2 * endPoint.y) + basePoint.y) / 3 };
		glm::vec2 base2 = { ((2 * basePoint.x) + endPoint.x) / 3 , ((2 * basePoint.y) + endPoint.y) / 3 };
		if (rec == 2) {
			std::vector<GLfloat> MahPoints = {
				startPoint.x, startPoint.y, -1,
				0,0,0,0,
				firstThird.x, firstThird.y, -1,
				0,0,0,0,
				midPoint.x, midPoint.y, -1,
				0,0,0,0,
				secondThird.x, secondThird.y, -1,
				0,0,0,0,
				endPoint.x, endPoint.y, -1,
				0,0,0,0
			};
			return MahPoints;
		}
		std::vector<GLfloat> segment1 = getMahPointsoutline(startPoint, firstThird, base2, rec - 1);
		std::vector<GLfloat> segment2 = getMahPointsoutline(firstThird, midPoint, secondThird, rec - 1);
		std::vector<GLfloat> segment3 = getMahPointsoutline(midPoint, secondThird, firstThird, rec - 1);
		std::vector<GLfloat> segment4 = getMahPointsoutline(secondThird, endPoint, base1, rec - 1);
		segment1.insert(segment1.end(), segment2.begin(), segment2.end());
		segment1.insert(segment1.end(), segment3.begin(), segment3.end());
		segment1.insert(segment1.end(), segment4.begin(), segment4.end());
		return segment1;
	}

	std::vector<GLfloat> snowflakeoutline(int rec) {
		glm::vec2 p0 = { -0.866 / 2.0f, -0.5f / 2.0f };
		glm::vec2 p1 = { 0, 0.5f };
		glm::vec2 p2 = { 0.866f / 2.0f, -0.5f / 2.0f };
		if (rec == 1) {
			std::vector<GLfloat> triangle = {
				p0.x,p0.y,-1,
				0,0,0,0,
				p1.x,p1.y,-1,
				0,0,0,0,
				p2.x,p2.y,-1,
				0,0,0,0
			};
			return triangle;
		}
		std::vector<GLfloat> points1 = getMahPointsoutline(p0, p1, p2, rec);
		std::vector<GLfloat> points2 = getMahPointsoutline(p1, p2, p0, rec);
		std::vector<GLfloat> points3 = getMahPointsoutline(p2, p0, p1, rec);
		points1.insert(points1.end(), points2.begin(), points2.end());
		points1.insert(points1.end(), points3.begin(), points3.end());
		return points1;
	}

	//gets all the points of the current recursive step, adds them to a vector and returns that vector

	std::vector<GLfloat> getMahPoints(glm::vec2 startPoint, glm::vec2 endPoint, glm::vec2 basePoint, int rec) {
		glm::vec2 firstThird = { ((2 * startPoint.x) + endPoint.x) / 3, ((2 * startPoint.y) + endPoint.y) / 3 };
		glm::vec2 secondThird = { ((2 * endPoint.x) + startPoint.x) / 3, ((2 * endPoint.y) + startPoint.y) / 3 };
		GLfloat Ux = endPoint.x - startPoint.x;
		GLfloat Uy = endPoint.y - startPoint.y;
		GLfloat Vx = endPoint.x - startPoint.x;
		GLfloat Vy = startPoint.y - endPoint.y;
		GLfloat Qx = startPoint.x + (1.0f / 2.0f)*Ux + (sqrt(3.0f) / 6.0f)*Vy;
		GLfloat Qy = startPoint.y + (1.0f / 2.0f)*Uy + (sqrt(3.0f) / 5.0f)*Vx;
		glm::vec2 midPoint = { Qx, Qy };
		glm::vec2 base1 = { ((2 * endPoint.x) + basePoint.x) / 3, ((2 * endPoint.y) + basePoint.y) / 3 };
		glm::vec2 base2 = { ((2 * basePoint.x) + endPoint.x) / 3 , ((2 * basePoint.y) + endPoint.y) / 3 };
		std::vector<GLfloat> MahPoints = {
			firstThird.x, firstThird.y, -1,
			1,0,0,1,
			midPoint.x, midPoint.y, -1,
			1,0,0,1,
			secondThird.x, secondThird.y, -1,
			1,0,0,1
		};
		if (rec == 2) {

			return MahPoints;
		}
		std::vector<GLfloat> segment1 = getMahPoints(startPoint, firstThird, base2, rec - 1);
		std::vector<GLfloat> segment2 = getMahPoints(firstThird, midPoint, secondThird, rec - 1);
		std::vector<GLfloat> segment3 = getMahPoints(midPoint, secondThird, firstThird, rec - 1);
		std::vector<GLfloat> segment4 = getMahPoints(secondThird, endPoint, base1, rec - 1);
		MahPoints.insert(MahPoints.end(), segment1.begin(), segment1.end());
		MahPoints.insert(MahPoints.end(), segment2.begin(), segment2.end());
		MahPoints.insert(MahPoints.end(), segment3.begin(), segment3.end());
		MahPoints.insert(MahPoints.end(), segment4.begin(), segment4.end());
		return MahPoints;
	}

	std::vector<GLfloat> snowflake(int rec) {
		glm::vec2 p0 = { -0.866 / 2.0f, -0.5f / 2.0f };
		glm::vec2 p1 = { 0, 0.5f };
		glm::vec2 p2 = { 0.866f / 2.0f, -0.5f / 2.0f };
		std::vector<GLfloat> triangle = {
			p0.x,p0.y,-1,
			1,0,0,1,
			p1.x,p1.y,-1,
			1,0,0,1,
			p2.x,p2.y,-1,
			1,0,0,1
		};
		if (rec == 1) {
			return triangle;
		}
		std::vector<GLfloat> points1 = getMahPoints(p0, p1, p2, rec);
		std::vector<GLfloat> points2 = getMahPoints(p1, p2, p0, rec);
		std::vector<GLfloat> points3 = getMahPoints(p2, p0, p1, rec);
		triangle.insert(triangle.end(), points1.begin(), points1.end());
		triangle.insert(triangle.end(), points2.begin(), points2.end());
		triangle.insert(triangle.end(), points3.begin(), points3.end());
		return triangle;
	}
	//void spinSnowflake(GLuint prog, GLfloat angle) {
	//	glUseProgram(prog);
	//	glm::mat4 trans = glm::mat4(1.0f);
	//	trans = glm::rotate(trans, angle, glm::vec3(0.0f, 0.0f, 1.0f));
	//	GLint uniTrans = glGetUniformLocation(prog, "rotationMtrx");
	//	glUniformMatrix4fv(uniTrans, 1, GL_FALSE, &trans[0][0]);
	//	glUseProgram(0);
	//}
	//-------------------------------------------- LAB 2 --------------------------------------------------------------------------------------------------------------


	struct triangle {
		glm::vec2 cI;
		glm::vec2 cJ;
		glm::vec2 c;
	};
	// Inte av mig. utan av sopan anton
	struct node {
		virtual ~node() { };
		node* parent = NULL;
	};

	struct BinaryNode :node {
		glm::vec2 c;
		glm::vec2 cI;
		glm::vec2 cM;
		glm::vec2 cJ;
		node* left = NULL;
		node* right = NULL;
	};

	struct ternaryNode :node {
		node* left = NULL;
		node* right = NULL;
		node* middle = NULL;
		glm::vec2 c;
		glm::vec2 cI;
		glm::vec2 cM;
		glm::vec2 cJ;
	};

	struct Leaf :node {
		triangle* tri;
	};

	std::vector<GLfloat> readFromfile() {
		std::vector<GLfloat> pointsVec;
		bool firstLine = true;
		//char c;
		int endOfLineCounter = 0;
		std::fstream file("points.txt", std::ios_base::in);
		float number;
		while (file >> number) {
			if (firstLine) {
				firstLine = false;
				continue;
			}
			else {
				endOfLineCounter++;
				pointsVec.push_back(number);
				if (endOfLineCounter == 2) {
					// z value
					pointsVec.push_back(-1);
					//colors
					pointsVec.push_back(1);
					pointsVec.push_back(0);
					pointsVec.push_back(0);
					pointsVec.push_back(1);
					endOfLineCounter = 0;
				}
			}

		}
		file.close();
		return pointsVec;
	}

	std::vector<GLfloat> randomizePoints(int numOfPoints) {
		std::vector<GLfloat> randomPoints;
		srand(static_cast <unsigned> (time(0)));
		for (int i = 0; i < numOfPoints; i++) {
			float x_coord = static_cast<float>(rand()) / (static_cast <float> (RAND_MAX / (1.9)));
			float x = x_coord - 0.9;
			randomPoints.push_back(x);
			GLfloat y_coord = static_cast<float>(rand()) / (static_cast <float> (RAND_MAX / (1.9)));
			float y = y_coord - 0.9;
			randomPoints.push_back(y);
			randomPoints.push_back(-1);  //z-value
			//color red
			randomPoints.push_back(1);
			randomPoints.push_back(0);
			randomPoints.push_back(0);
			randomPoints.push_back(1);
		}
		return randomPoints;
	}

	std::vector<glm::vec2> coordinates(std::vector<GLfloat> pointVector) {
		std::vector<glm::vec2> coordVec;
		int size = pointVector.size() / 7;
		int it = 0;
		//std::cout << size << std::endl;
		for (int i = 0; i < pointVector.size(); i = i + 7) {
			it = it + 1;
			glm::vec2 point(pointVector[i], pointVector[i + 1]);
			coordVec.push_back(point);
			if (it == size) {
				//	std::cout << coordVec[1][1] << std::endl;
					//printf("Hej");
				return coordVec;
			}
		}
	}

	std::vector<GLfloat> vec3ToGLfloat(std::vector<glm::vec2> vector) {
		std::vector<GLfloat> vectorGLfloat;
		for (int i = 0; i < vector.size(); i++) {
			vectorGLfloat.push_back(vector[i][0]);
			vectorGLfloat.push_back(vector[i][1]);
			vectorGLfloat.push_back(-1);
			vectorGLfloat.push_back(1);
			vectorGLfloat.push_back(0);
			vectorGLfloat.push_back(0);
			vectorGLfloat.push_back(1);
			//	std::cout << vector[i][0] << std::endl;
			//	std::cout << vector[i][1] << std::endl;
			//	std::cout << vector[i][2] << std::endl;
		}
		return vectorGLfloat;
	}

	inline void sortPoints(std::vector<glm::vec2>& points) {
		std::sort(
			points.begin(),
			points.end(),
			[](const glm::vec2& v1, const glm::vec2& v2) {
			return (v1.x == v2.x) ? (v1.y < v2.y) : (v1.x < v2.x);
		}
		);
	}


	bool checkClockwise1(glm::vec2 A, glm::vec2 B, glm::vec2 C) {
		return ((B.x - A.x) * (C.y - A.y)) > ((B.y - A.y) * (C.x - A.x));
	}

	int turnClockwise(glm::vec2 startPoint, glm::vec2 endPoint, glm::vec2 newPoint) {
		float num = ((newPoint.x - startPoint.x) * (endPoint.y - startPoint.y)) - ((newPoint.y - startPoint.y)*(endPoint.x - startPoint.x));
		if (num < 0) {
			return 1;
		}
		//else if (num == 0) {
		//	return 0;
		//}
		else {
			return -1;
		}
	}

	std::vector<glm::vec2> andrews(std::vector<GLfloat> setOfPoints) {
		std::vector<glm::vec2> points = coordinates(setOfPoints);
		sortPoints(points);
		if (points.size() <= 3) {
			return points;
		}
		//	std::cout << points.size() << std::endl;
		std::vector<glm::vec2> upper;
		std::vector<glm::vec2> lower;
		std::vector<glm::vec2> hull;
		for (int i = 0; i < points.size(); i++) {
			int sizeL = lower.size();
			while (sizeL >= 2 && checkClockwise1(points.at(i), lower.at(sizeL - 1), lower.at(sizeL - 2))) {
				lower.pop_back();
				sizeL = lower.size();
			}
			lower.push_back(points.at(i));
		}
		for (int i = points.size() - 1; i >= 0; i--) {
			int sizeU = upper.size();
			while (sizeU >= 2 && checkClockwise1(points.at(i), upper.at(sizeU - 1), upper.at(sizeU - 2))) {
				upper.pop_back();
				sizeU = upper.size();
			}
			upper.push_back(points[i]);
		}
		upper.pop_back();
		hull.insert(hull.end(), upper.begin(), upper.end());
		hull.insert(hull.end(), lower.begin(), lower.end());
		return hull;
	}


	glm::vec2 findMedianVertex(std::vector<glm::vec2> vector) {
		int sizeOfVec = vector.size();
		return vector.at(sizeOfVec / 2);
	}

	std::vector<glm::vec2> removeHull(std::vector<glm::vec2> pointSet, std::vector<glm::vec2> convexHull) {
		for (int i = 0; i < convexHull.size(); i++) {
			for (int j = 0; j < pointSet.size(); j++) {
				if (convexHull.at(i) == pointSet.at(j)) {
					pointSet.erase(pointSet.begin() + j);
				}
				else {
					continue;
				}
			}
		}
		return pointSet;
	}

	glm::vec2 selectRandomPoint(std::vector<glm::vec2> points) {
		return points.at(points.size() / 2);
	}

	node* buildTree(glm::vec2 c, std::vector<glm::vec2> convexHull, node* parent) {
		glm::vec2 median = findMedianVertex(convexHull);
		if (convexHull.size() == 2 && parent != nullptr) {
			Leaf* leafNode = new Leaf();
			triangle* trngl = new triangle();
			trngl->c = c;
			trngl->cI = convexHull.at(0);
			trngl->cJ = convexHull.at(1);
			leafNode->tri = trngl;
			leafNode->parent = parent;
			return leafNode;
		}
		BinaryNode* b = new BinaryNode();
		int medianLocation = convexHull.size() / 2;
		std::vector<glm::vec2> firstHalf(convexHull.begin(), convexHull.begin() + medianLocation + 1);
		//std::cout << firstHalf.size() << std::endl;
		std::vector<glm::vec2> secondHalf(convexHull.begin() + medianLocation, convexHull.end());
		node* L = new node();
		node* R = new node();
		L = buildTree(c, firstHalf, b);
		R = buildTree(c, secondHalf, b);
		glm::vec2 cI = convexHull.front();
		glm::vec2 cJ = convexHull.back();
		b->left = L;
		b->right = R;
		b->parent = parent;
		b->c = c;
		b->cI = cI;
		b->cJ = cJ;
		b->cM = median;
		return b;
	}

	std::vector<glm::vec2> treeToList(node* tree) {
		std::vector<glm::vec2> lineList;
		if (Leaf* leaf = dynamic_cast<Leaf*>(tree)) {
			triangle* trngl = leaf->tri;
			lineList.push_back(trngl->c);
			lineList.push_back(trngl->cI);

			lineList.push_back(trngl->cI);
			lineList.push_back(trngl->cJ);

			lineList.push_back(trngl->cJ);
			lineList.push_back(trngl->c);
			return lineList;
		}
		else if (BinaryNode* binarynode = dynamic_cast<BinaryNode*>(tree)) {
			std::vector<glm::vec2> leftList = treeToList(binarynode->left);
			std::vector<glm::vec2> rightList = treeToList(binarynode->right);
			lineList.insert(lineList.end(), leftList.begin(), leftList.end());
			lineList.insert(lineList.end(), rightList.begin(), rightList.end());
			return lineList;
		}
		else if (ternaryNode* ternary = dynamic_cast<ternaryNode*>(tree)) {
			std::vector<glm::vec2> leftList = treeToList(ternary->left);
			std::vector<glm::vec2> middleList = treeToList(ternary->middle);
			std::vector<glm::vec2> rightList = treeToList(ternary->right);
			lineList.insert(lineList.end(), leftList.begin(), leftList.end());
			lineList.insert(lineList.end(), middleList.begin(), middleList.end());
			lineList.insert(lineList.end(), rightList.begin(), rightList.end());
			return lineList;
		}
	}

	float checkSector(glm::vec2 start, glm::vec2 end, glm::vec2 point) {
		return ((end.x - start.x) * (point.y - start.y)) - ((end.y - start.y) * (point.x - start.x));
	}

	int pointLocation(glm::vec2 leftPoint, glm::vec2 midPoint, glm::vec2 rightPoint, glm::vec2 q) {
		int thecase = turnClockwise(midPoint, leftPoint, rightPoint);
		int rightOfcic = turnClockwise(midPoint, rightPoint, q);
		int leftOfcmc = turnClockwise(midPoint, leftPoint, q);
		if (thecase < 1) {
			if (rightOfcic == 1 && leftOfcmc == -1) {
				return 1;
			}
			//else if ((rightOfcic == 0 && leftOfcmc == 0) || (rightOfcic == 0 && leftOfcmc == -1) || (rightOfcic == 1 && leftOfcmc == 0)) {
			//	return 0;
			//}
			else {
				return -1;
			}

		}

		else {
			if (rightOfcic == 1 || leftOfcmc == -1) {
				return 1;
			}
			//else if ((rightOfcic == 0 || leftOfcmc == 0) || (rightOfcic == 0 || leftOfcmc == -1) || (rightOfcic == 1 || leftOfcmc == 0)) {
			//	return 0;
			//}
			else {
				return -1;
			}
		}

	}

	node* insertOnTheLine(node* tree, glm::vec2 onLinePoint) {
		if (Leaf* leaf = dynamic_cast<Leaf*>(tree)) {
			BinaryNode* binNode = new BinaryNode();

			Leaf* leafRight = new Leaf();
			Leaf* leafLeft = new Leaf();

			triangle* triLeft = new triangle();
			triangle* triRight = new triangle();

			glm::vec2 c = leaf->tri->c;
			glm::vec2 ci = leaf->tri->cI;
			glm::vec2 cj = leaf->tri->cJ;

			int isOnLineCJCI = turnClockwise(cj, ci, onLinePoint);
			int isOnLineCJC = turnClockwise(cj, c, onLinePoint);
			int isOnLineCCI = turnClockwise(c, ci, onLinePoint);
			if (isOnLineCJCI == 0) {
				triLeft->c = onLinePoint;
				triLeft->cI = c;
				triLeft->cJ = cj;

				triRight->c = onLinePoint;
				triRight->cI = ci;
				triRight->cJ = c;
			}
			else if (isOnLineCCI == 0) {
				triLeft->c = onLinePoint;
				triLeft->cI = ci;
				triLeft->cJ = cj;

				triRight->c = onLinePoint;
				triRight->cI = cj;
				triRight->cJ = c;
			}
			else {
				triLeft->c = onLinePoint;
				triLeft->cI = cj;
				triLeft->cJ = ci;

				triRight->c = onLinePoint;
				triRight->cI = ci;
				triRight->cJ = c;
			}
			leafLeft->tri = triLeft;
			leafRight->tri = triRight;

			binNode->cI = ci;
			binNode->cJ = cj;
			binNode->cM = c;
			binNode->c = onLinePoint;
			binNode->left = leafLeft;
			binNode->right = leafRight;
			binNode->parent = leaf->parent;
			return binNode;
		}
		else if (BinaryNode* binNode = dynamic_cast<BinaryNode*>(tree)) {
			glm::vec2 c = binNode->c;
			glm::vec2 ci = binNode->cI;
			glm::vec2 cj = binNode->cJ;
			glm::vec2 cm = binNode->cM;
			int searchLeft = pointLocation(cm, c, ci, onLinePoint);
			if (searchLeft == 0) {
				binNode->left = insertOnTheLine(binNode->left, onLinePoint);
				binNode->left->parent = binNode;
				return binNode;
			}
			else {
				binNode->right = insertOnTheLine(binNode->right, onLinePoint);
				binNode->right->parent = binNode;
				return binNode;
			}
		}
		else if (ternaryNode* ternary = dynamic_cast<ternaryNode*>(tree)) {
			glm::vec2 c = ternary->c;
			glm::vec2 ci = ternary->cI;
			glm::vec2 cj = ternary->cJ;
			glm::vec2 cm = ternary->cM;

			int searchLeftTernary = pointLocation(ci, c, cm, onLinePoint);
			int searchMidTernary = pointLocation(cj, c, ci, onLinePoint);

			int inLeftBase = turnClockwise(ci, cj, onLinePoint);
			int inMidBase = turnClockwise(ci, cm, onLinePoint);
			int inRightBase = turnClockwise(cj, cm, onLinePoint);
			if (searchLeftTernary == 0 || inLeftBase == 0) { //&& eller || ???
				ternary->left = insertOnTheLine(ternary->left, onLinePoint);
				ternary->left->parent = ternary;
				return ternary;
			}
			else if (searchMidTernary == 0 || inMidBase == 0) {
				ternary->middle = insertOnTheLine(ternary->middle, onLinePoint);
				ternary->middle->parent = ternary;
				return ternary;
			}
			else {
				ternary->right = insertOnTheLine(ternary->right, onLinePoint);
				ternary->right->parent = ternary;
				return ternary;
			}
		}
	}

	node* insert(node* tree, glm::vec2 point) {
		if (Leaf* leaf = dynamic_cast<Leaf*>(tree)) {
			triangle* leafTri = leaf->tri;

			triangle* triLeft = new triangle();
			triangle* triMid = new triangle();
			triangle* triRight = new triangle();
			Leaf* leafLeft = new Leaf();
			Leaf* leafMid = new Leaf();
			Leaf* leafRight = new Leaf();

			glm::vec2 ci = leafTri->cI;
			glm::vec2 cj = leafTri->cJ;
			glm::vec2 c = leafTri->c;

			triLeft->cI = c;
			triLeft->c = point;
			triLeft->cJ = ci;
			leafLeft->tri = triLeft;

			triMid->cI = ci;
			triMid->c = point;
			triMid->cJ = cj;
			leafMid->tri = triMid;

			triRight->cI = cj;
			triRight->c = point;
			triRight->cJ = c;
			leafRight->tri = triRight;

			ternaryNode* ternary = new ternaryNode();

			leafRight->parent = ternary;
			leafMid->parent = ternary;
			leafLeft->parent = ternary;

			ternary->left = leafLeft;
			ternary->middle = leafMid;
			ternary->right = leafRight;

			ternary->c = point;
			ternary->cI = ci;
			ternary->cJ = cj;
			ternary->cM = c;

			ternary->parent = leaf->parent;
			return ternary;
		}
		else if (BinaryNode* binNode = dynamic_cast<BinaryNode*>(tree)) {
			glm::vec2 c = binNode->c;
			glm::vec2 ci = binNode->cI;
			glm::vec2 cm = binNode->cM;
			glm::vec2 cj = binNode->cJ;
			int locationOfPoint = pointLocation(cm, c, ci, point);
			int loctionOfPointRight = pointLocation(cj, c, cm, point);
			if (locationOfPoint == 1) {
				binNode->left = insert(binNode->left, point);
				binNode->left->parent = binNode;
				return binNode;
			}
			else if (locationOfPoint == -1) {
				binNode->right = insert(binNode->right, point);
				binNode->right->parent = binNode;
				return binNode;
			}
		}
		else if (ternaryNode* ternary = dynamic_cast<ternaryNode*>(tree)) {
			glm::vec2 c = ternary->c;
			glm::vec2 ci = ternary->cI;
			glm::vec2 cm = ternary->cM;
			glm::vec2 cj = ternary->cJ;
			int locationOfPointLeft = pointLocation(ci, c, cm, point);
			int locationOfPointMiddle = pointLocation(cj, c, ci, point);
			int locationOfPointRight = pointLocation(cm, c, cj, point);
			if (locationOfPointLeft == 1) {
				ternary->left = insert(ternary->left, point);
				ternary->left->parent = ternary;
				return ternary;
			}
			else if (locationOfPointLeft == 0) {
				ternary->left = insertOnTheLine(ternary->left, point);
				ternary->left->parent = ternary;
				return ternary;
			}
			else if (locationOfPointMiddle == 1) {
				ternary->middle = insert(ternary->middle, point);
				ternary->middle->parent = ternary;
				return ternary;
			}
			else if (locationOfPointMiddle == 0) {
				ternary->middle = insertOnTheLine(ternary->middle, point);
				ternary->middle->parent = ternary;
				return ternary;
			}
			else if (locationOfPointRight == 1) {
				ternary->right = insert(ternary->right, point);
				ternary->right->parent = ternary;
				return ternary;
			}
			else {
				ternary->right = insertOnTheLine(ternary->right, point);
				ternary->right->parent = ternary;
				return ternary;
			}
		}
	}

	node* getTriangulationSoup(node* tree, std::vector<glm::vec2> pointsInsideHull) {
		for (int i = 0; i < pointsInsideHull.size(); i++) {
			glm::vec2 p = pointsInsideHull.at(i);
			tree = insert(tree, p);
		}
		return tree;
	}

	size_t bufSize;
	size_t bufSizeAndrews;
	size_t bufSizeTriangleFan;
	size_t bufSizeTriangulation;
	bool

		ExampleApp::Open()
	{
		App::Open();



		this->window = new Display::Window;
		window->SetKeyPressFunction([this](int32 key, int32, int32 action, int32)
		{
			if ((GLFW_KEY_UP == key) && (action == GLFW_RELEASE)) {
				randomizePoints(8);
			}
			else if ((GLFW_KEY_DOWN == key) && (action == GLFW_RELEASE)) {
				if (rec > 1) {
					rec--;
				}

			}
			else if ((GLFW_KEY_ESCAPE == key) && (action == GLFW_RELEASE) || (rec < 1)) {
				this->window->Close();
			}
		});

		//std::vector<GLfloat> baseTri = readFromfile();
		//std::cout << baseTri.size() << std::endl;

		//bufSize = baseTri.size() / 7;
		//GLfloat* buf = &baseTri[0];


		if (this->window->Open())
		{
			// set clear color to gray
			glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
			this->vsBuffer = new GLchar[STRING_BUFFER_SIZE];
			this->fsBuffer = new GLchar[STRING_BUFFER_SIZE];

			// copy the hardcoded shader into buffer
			strncpy_s(this->vsBuffer, STRING_BUFFER_SIZE, vs, STRING_BUFFER_SIZE);
			strncpy_s(this->fsBuffer, STRING_BUFFER_SIZE, ps, STRING_BUFFER_SIZE);

			// setup vertex shader
			this->vertexShader = glCreateShader(GL_VERTEX_SHADER);
			GLint length = (GLint)std::strlen(vs);
			glShaderSource(this->vertexShader, 1, &vs, &length);
			glCompileShader(this->vertexShader);

			// get error log
			GLint shaderLogSize;
			glGetShaderiv(this->vertexShader, GL_INFO_LOG_LENGTH, &shaderLogSize);
			if (shaderLogSize > 0)
			{
				GLchar* buf = new GLchar[shaderLogSize];
				glGetShaderInfoLog(this->vertexShader, shaderLogSize, NULL, buf);
				printf("[SHADER COMPILE ERROR]: %s", buf);
				delete[] buf;
			}

			// setup pixel shader
			this->pixelShader = glCreateShader(GL_FRAGMENT_SHADER);
			length = (GLint)std::strlen(ps);
			glShaderSource(this->pixelShader, 1, &ps, &length);
			glCompileShader(this->pixelShader);

			// get error log
			shaderLogSize;
			glGetShaderiv(this->pixelShader, GL_INFO_LOG_LENGTH, &shaderLogSize);
			if (shaderLogSize > 0)
			{
				GLchar* buf = new GLchar[shaderLogSize];
				glGetShaderInfoLog(this->pixelShader, shaderLogSize, NULL, buf);
				printf("[SHADER COMPILE ERROR]: %s", buf);
				delete[] buf;
			}

			// create a program object
			this->program = glCreateProgram();
			glAttachShader(this->program, this->vertexShader);
			glAttachShader(this->program, this->pixelShader);
			glLinkProgram(this->program);
			glGetProgramiv(this->program, GL_INFO_LOG_LENGTH, &shaderLogSize);
			if (shaderLogSize > 0)
			{
				GLchar* buf = new GLchar[shaderLogSize];
				glGetProgramInfoLog(this->program, shaderLogSize, NULL, buf);
				printf("[PROGRAM LINK ERROR]: %s", buf);
				delete[] buf;
			}

			// setup vbo
			glGenBuffers(1, &this->triangle);
			glBindBuffer(GL_ARRAY_BUFFER, this->triangle);
			glEnableVertexAttribArray(0);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float32) * 7, NULL);
			glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(float32) * 7, (GLvoid*)(sizeof(float32) * 3));
			//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*baseTri.size(), baseTri.data(), GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			this->window->SetUiRender([this]()
			{
				this->RenderUI();
			});
			return true;
		}
		return false;
	}

	//------------------------------------------------------------------------------
	/**
	*/
	//GLfloat angle = 0;
	void
		ExampleApp::Run()
	{
		glLineWidth(3.0f);
		glPointSize(10.0f);
		std::vector<GLfloat> baseTri = readFromfile();
		std::vector<glm::vec2> triPoints = coordinates(baseTri);
		bufSize = baseTri.size() / 7;

		std::vector<glm::vec2> hullPoints = andrews(baseTri);
		std::vector<glm::vec2> insideHullPoints = removeHull(triPoints, hullPoints);

		std::vector<GLfloat> hull = vec3ToGLfloat(hullPoints);
		for (int i = 0; i < hull.size(); i++) {
			std::cout << hull.at(i) << std::endl;
		}

		node* startNode = NULL;
		startNode = new node();
		node* tree = new node();
		glm::vec2 randomPoint = selectRandomPoint(insideHullPoints);
		tree = buildTree(randomPoint, hullPoints, startNode);
		std::vector<glm::vec2> treeList = treeToList(tree);
		std::vector<GLfloat> floatPointList = vec3ToGLfloat(treeList);
		bufSizeTriangleFan = floatPointList.size() / 7;
		//std::cout << floatPointList.size() << std::endl;

		node* fullTree = new node();
		fullTree = getTriangulationSoup(tree, insideHullPoints);
		std::vector<glm::vec2> extendedTree = treeToList(fullTree);
		std::vector<GLfloat> extendedTreeGLfloat = vec3ToGLfloat(extendedTree);
		bufSizeTriangulation = extendedTreeGLfloat.size() / 7;
		//std::cout << extendedTree.size() << std::endl;
		//std::cout << extendedTree.size() << std::endl;
		bufSizeAndrews = hull.size();
		while (this->window->IsOpen())
		{


			glClear(GL_COLOR_BUFFER_BIT);
			this->window->Update();
			//spinSnowflake(this->program, glm::radians(angle));
			//angle -= 1;
			// setup vbo
			glUseProgram(this->program);
			glBindBuffer(GL_ARRAY_BUFFER, this->triangle);
			glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*extendedTreeGLfloat.size(), extendedTreeGLfloat.data(), GL_STATIC_DRAW);
			glDrawArrays(GL_LINES, 0, bufSizeTriangulation);
			glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*floatPointList.size(), floatPointList.data(), GL_STATIC_DRAW);
			glDrawArrays(GL_LINES, 0, bufSizeTriangleFan);
			glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*hull.size(), hull.data(), GL_STATIC_DRAW);
			glDrawArrays(GL_LINE_LOOP, 0, bufSizeAndrews / 7);
			glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*baseTri.size(), baseTri.data(), GL_STATIC_DRAW);
			glDrawArrays(GL_POINTS, 0, bufSize);
			glUseProgram(0);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			// Outline -----------------------------
	//		std::vector<glm::vec3> baseTri1 = andrews(baseTri);

	//		bufSize = baseTri1.size() / 7;
			//buf = &baseTri[0];

			//// setup vbo
			//glGenBuffers(1, &this->triangle);
			//glBindBuffer(GL_ARRAY_BUFFER, this->triangle);
			//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*baseTri1.size(), baseTri1.data(), GL_STATIC_DRAW);
			//glBindBuffer(GL_ARRAY_BUFFER, 0);

			// do stuff
			//glBindBuffer(GL_ARRAY_BUFFER, this->triangle);
			//glUseProgram(this->program);
			//glEnableVertexAttribArray(0);
			//glEnableVertexAttribArray(1);
			//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float32) * 7, NULL);
			//glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(float32) * 7, (GLvoid*)(sizeof(float32) * 3));
			//glDrawArrays(GL_LINE_LOOP, 0, bufSize);
			//glBindBuffer(GL_ARRAY_BUFFER, 0);

			this->window->SwapBuffers();
		}
	}
	void
		ExampleApp::RenderUI()
	{
		if (this->window->IsOpen())
		{
			bool show = true;
			static int sliderInt = 6;
			// create a new window
			ImGui::Begin("Shader Sources", &show, ImGuiWindowFlags_NoSavedSettings);

			// create text editors for shader code

			// apply button
			if (ImGui::Button("Read from file"))
			{
			}
			else if (ImGui::Button("Randomize")) {

			}
			ImGui::SliderInt("amount of points", &sliderInt, 6, 50);
			if (this->compilerLog.length())
			{
				// if compilation produced any output we display it here
				ImGui::TextWrapped(this->compilerLog.c_str());
			}
			// close window
			ImGui::End();
		}
	}
} // namespace Example

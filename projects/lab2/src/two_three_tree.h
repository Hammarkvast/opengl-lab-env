#include <glm/vec2.hpp>
#include <iostream>
#include "types.h"
#include "point_line.h"

#ifndef TWO_THREE_TREE_H
#define TWO_THREE_TREE_H

struct Triangle {
	Point p0, p1, p2;

	Triangle(Point p0, Point p1, Point p2) :
		p0(p0), p1(p1), p2(p2) {}

	PointSet toVec() {
		return {p0, p1, p2};
	}
};

struct Node {
	Node *parent;

	virtual ~Node(){}
	virtual Point getC(){}
	virtual Point getCm(){}
	virtual Point getCi(){}
	virtual Point getCj(){}

	virtual Node* getLst(){}
	virtual Node* getMst(){}
	virtual Node* getRst(){}
	virtual void setLst(Node* tree){}
	virtual void setMst(Node* tree){}
	virtual void setRst(Node* tree){}

	virtual void insertPoint(Point &point) = 0;

	Node(Node *parent) :
		parent(parent) {};
};

struct BNode : Node {
	Point c, ci, cm, cj;
	Node *lst, *rst;

	BNode(Point c, Point ci, Point cm, Point cj, Node *parent) :
		c(c), ci(ci), cm(cm), cj(cj), lst(nullptr), rst(nullptr), Node(parent) {}
	
	virtual Point getC() {return c;}
	virtual Point getCm() {return cm;}
	virtual Point getCi() {return ci;}
	virtual Point getCj() {return cj;}

	virtual Node* getLst() {return lst;}
	virtual Node* getMst() {return nullptr;}
	virtual Node* getRst() {return rst;}
	virtual void setLst(Node* tree) {lst = tree;}
	virtual void setRst(Node* tree) {rst = tree;}

	PointSet toVec() {
		return {c, ci, cm, cj}; 
	}

	void insertPoint(Point &point) {
		printf("SEARCHING IN BINARY NODE\n");
		printf("point	: (%f, %f)\n", point.x, point.y);
		printf("c	: (%f, %f)\n", c.x, c.y);
		printf("ci	: (%f, %f)\n", ci.x, ci.y);
		printf("cm	: (%f, %f)\n", cm.x, cm.y);
		printf("cj	: (%f, %f)\n", cj.x, cj.y);

		if (leftOf(cm, c, ci) || onLine(cm, c, ci, false)) {	// case 1 "ci left of cm->c"
			if (!leftOf(ci, c, point) && leftOf(cm, c, point)) {
				printf("left case1\n");
				lst->insertPoint(point);
			} else {
				printf("in right case1\n");
				rst->insertPoint(point);
			}
		} else { 	// case 2 "ci right of cm->c"
			if (!leftOf(ci, c, point) || leftOf(cm, c, point)) {
				printf("left case2\n");
				lst->insertPoint(point);
			} else {
				printf("right case2\n");
				rst->insertPoint(point);
			}
		}
	}
};

struct TNode : Node {
	Point c, ci, cm, cj;
	Node *lst, *mst, *rst;

	TNode(Point c, Point ci, Point cm, Point cj, Node *parent) :
		c(c), ci(ci), cm(cm), cj(cj), lst(nullptr), mst(nullptr), rst(nullptr), Node(parent) {}
	
	virtual Point getC() {return c;}
	virtual Point getCm() {return cm;}
	virtual Point getCi() {return ci;}
	virtual Point getCj() {return cj;}

	virtual Node* getLst() {return lst;}
	virtual Node* getMst() {return mst;}
	virtual Node* getRst() {return rst;}
	virtual void setLst(Node* tree) {lst = tree;}
	virtual void setMst(Node* tree) {mst = tree;}
	virtual void setRst(Node* tree) {rst = tree;}

	virtual void insertPoint(Point &point) {
		return;
	}

	PointSet toVec() {
		return {c, ci, cm, cj};
	}
};

struct Leaf : Node {
	Triangle *triangle;

	Leaf(Triangle *triangle, Node *parent) :
		triangle(triangle), Node(parent) {}

	void insertPoint(Point &point) {
		TNode *newNode = new TNode(point, triangle->p0, triangle->p1, triangle->p2, parent);
		newNode->lst = new Leaf(new Triangle(triangle->p0, triangle->p1, point), newNode);
		newNode->mst = new Leaf(new Triangle(point, triangle->p1, triangle->p2), newNode);
		newNode->rst = new Leaf(new Triangle(triangle->p0, point, triangle->p2), newNode);

		printf("FOUND IN LEAF\n");
		printf("(%f, %f) -> (%f, %f) -> (%f, %f)\n", triangle->p0.x,triangle->p0.y,triangle->p1.x,triangle->p1.y,triangle->p2.x,triangle->p2.y);

		// Determine the sub-tree of parent leaf is located
		if (parent->getLst() == this) {
			parent->setLst(newNode);
		} else if (parent->getMst() == this) {
			parent->setMst(newNode);
		} else {
			parent->setRst(newNode);
		}
	}
};

#endif
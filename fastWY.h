#ifndef FASTWY_H
#define FASTWY_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <boost/math/distributions/hypergeometric.hpp>
#include <set>

#include "database.h"

using namespace std;
using namespace boost::math;

using uint = unsigned int;
using ullint = unsigned long long int;
using uchar = unsigned char;

class FastWY{
private:

	struct Position{
		uint sequence_index;
		uint event_index;
		uint itemset_index;
		bool operator==(const Position& x){
			return (sequence_index == x.sequence_index && event_index == x.event_index && itemset_index == x.itemset_index);
		}
		bool operator!=(const Position& x){
			return (sequence_index != x.sequence_index || event_index != x.event_index || itemset_index != x.itemset_index);
		}
	};

	using projectDB = vector< Position >;

	struct Node{
		list<Node>::iterator parent;
		vector<list<Node>::iterator> child;
		vector<Event> pattern;
		string patternStr;
		uint supportSum;
		vector<uint> support;
		projectDB pdb;
		vector<uint> x;
		double f_val;
		double p;
		bool closed;
		uint pdb_size;
	};

	const uint MAXVAL = 0xffffffff;
	uint mItemSize;
	double mMinsup;
	uint mMinpat;
	uint mMaxpat;
	uint mMinItem;
	uint mMaxItem;
	int mInterval;
	uint mSupMode;
	uint mCloSpan;
	uint mWild;
	uint mSide;

	uint mN;
	vector<double> mY;
	double mDelta;

	vector<double> mL;
	uint mR;

	vector<vector<double>> mTempY;
	vector<double> mPmin;
	uint mPlusN;
	double mAlpha;

	using Key = pair<uint, uint>;
	map<Key, double> mPvalueMap;
	vector<vector<Event> > mTransaction;
	vector<Event> mPattern;
	list<Node> mTree;
	using TreeIter = list<Node>::iterator;
	unordered_map<uint, list<TreeIter>> mSignificant;

	uint mFlagCScheckEnd = 0;
	unordered_map<uint, list<TreeIter>> mCheckPDB;

	uint mFlagItemsetExist = 0;
	Event mWildEvent;

	string pat2str(void);
	bool calculate(Node& aNode);
	bool calculate_new(Node& aNode);
	void project(const TreeIter aNodeIter);
	void project_new(projectDB& aPdb, const TreeIter aParent);
	uint calcSup(uint aId, vector<Event> aPattern);
	uint calcPat(uint aId, vector<Event> aPattern);
	uint isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2);
	bool isInclude(const Event aEvent1, const Event aEvent2);
	bool isParent(const TreeIter aNodeIter, const TreeIter aChildIter);
	void checkProjectedDB(const TreeIter aCurrentIter);
	uint sum_items_pdb(const projectDB aPdb);
	uint sumId(const vector<uint> aIdList);
	double exactTest(const Node aNode);
	vector<double> exactTest_vector(const Node aNode);
	void childPatternUpdate(const TreeIter aChildIter);
	void search_tree(const TreeIter aNodeIter);
	void search_tree_final(const TreeIter aNodeIter);
	double log_combination(const uint aS, const uint aT);

public:
	FastWY(double aMinsup, uint aMinpat, uint aMaxpat, uint aMinItem, uint aMaxItem, uint aR, double aAlpha, uint aCloSpan, uint aSupMode, int aInterval, uint aWild, uint aSide){
		mMinsup = aMinsup;
		mMinpat = aMinpat;
		mMinItem = aMinItem;
		mMaxpat = aMaxpat;
		mMaxItem = aMaxItem;
		mR = aR;
		mAlpha = aAlpha;
		mCloSpan = aCloSpan;
		mSupMode = aSupMode;
		mInterval = aInterval;
		mWild = aWild;
		mSide = aSide;

	}
	;

	void init(const vector<vector<Event> > aTransaction, const vector<double> aY);
	void main();
	void printTree(string aFilename);
	void printSigPattern(string aFilename);
};

#endif

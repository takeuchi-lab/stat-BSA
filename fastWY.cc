#include "fastWY.h"


string FastWY::pat2str(void){
	stringstream ss;
	for(uint i = 0; i < mPattern.size(); ++i){
		ss << (i ? " " : "");
		if(mPattern[i].first.size() > 0){
			ss << "(";
			for(uint j = 0; j < mPattern[i].first.size(); ++j){
				ss << mPattern[i].first[j];
				ss << (j + 1 < mPattern[i].first.size() ? "_" : "");
			}
			ss << ")";
			if(mPattern[i].second.size() > 0) ss << ":";
		}
		for(uint j = 0; j < mItemSize; j++){
			if(mPattern[i].second[j] == MAXVAL) ss << "*";
			else ss << mPattern[i].second[j];
			ss << (j + 1 < mItemSize ? ":" : "");
		}
	}
	return ss.str();
}



void FastWY::printTree(string aFilename){
	ofstream tFile(aFilename);
	tFile << "pattern(d = " << mDelta << "),length,p,f,supportSum,sup+,sup-,+ or -" << '\n';
	for(auto it = ++mTree.begin(); it != mTree.end(); ++it){
		tFile << it->patternStr << "," << it->pattern.size() << "," << it->p << "," << it->f_val << "," << it->supportSum << ",";
		uint sup_p = 0;
		uint sup_m = 0;
		for(uint j = 0; j < it->support.size(); j++){
			if(mY[it->x[j]] > 0) sup_p += it->support[j];
			else sup_m += it->support[j];
		}
		tFile << sup_p << "," << sup_m << ",";
		if(sup_p > (double) it->supportSum * mPlusN / mN) tFile << "+";
		else if(sup_p < (double) it->supportSum * mPlusN / mN) tFile << "-";
		else tFile << "N";
		tFile << '\n';
	}
}


void FastWY::printSigPattern(string aFilename){
	ofstream tFile(aFilename);
	tFile << "pattern(d = " << mDelta << "),length,p,f,supportSum,sup+,sup-,+ or -" << "\n";
	for(auto it : mSignificant){
		for(auto itr : it.second){
			tFile << itr->patternStr << "," << itr->pattern.size() << "," << itr->p << "," << itr->f_val << "," << itr->supportSum << ",";
			uint sup_p = 0;
			uint sup_m = 0;
			for(uint j = 0; j < itr->support.size(); j++){
				if(mY[itr->x[j]] > 0) sup_p += itr->support[j];
				else sup_m += itr->support[j];
			}
			tFile << sup_p << "," << sup_m << ",";
			if(sup_p > (double) itr->supportSum * mPlusN / mN) tFile << "+";
			else if(sup_p < (double) itr->supportSum * mPlusN / mN) tFile << "-";
			else tFile << "N";

			tFile << "\n";
		}
	}
}


double FastWY::log_combination(const uint aS, const uint aT){
	double tLog_s = 0;
	double tLog_t = 0;

	for(uint i = aS; i > aS - aT; --i){
		tLog_s += log(i);
	}

	for(uint i = aT; i > 0; --i){
		tLog_t += log(i);
	}

	return tLog_s - tLog_t;
}



double FastWY::exactTest(const Node aNode){

	vector<uint> tX = aNode.x;

	uint count = 0;
	for(auto it : tX){
		if(mY[it] == 1) count++;
	}

	Key tKey = pair<uint, uint>(count, aNode.supportSum);
	if(mPvalueMap.count(tKey) == 0){
		uint max_for_k = (mSide != 2) ? min(mPlusN, aNode.supportSum) : count;
		uint min_for_k = (mSide != 1) ? (uint) max(0, int(mPlusN + aNode.supportSum - mN)) : count;
		hypergeometric_distribution<> hgd(mPlusN, aNode.supportSum, mN);
		double cutoff = pdf(hgd, count);
		double tmp_p = 0.0;

		for(uint k = min_for_k; k <= max_for_k; k++){
			double p = pdf(hgd, k);
			if(mSide == 0){
				if(p <= cutoff){
					tmp_p += p;
				}
			}else{
				tmp_p += p;
			}

		}
		mPvalueMap[tKey] = tmp_p;

	}
	return mPvalueMap[tKey];
}


vector<double> FastWY::exactTest_vector(const Node aNode){
	vector<double> tP;
	vector<uint> tX = aNode.x;
	for(uint i = 0; i < mR; ++i){

		uint count = 0;
		for(auto it : tX){
			if(mTempY[i][it] == 1) count++;
		}

		Key tKey = pair<uint, uint>(count, aNode.supportSum);
		if(mPvalueMap.count(tKey) == 0){
			uint max_for_k = (mSide != 2) ? min(mPlusN, aNode.supportSum) : count;
			uint min_for_k = (mSide != 1) ? (uint) max(0, int(mPlusN + aNode.supportSum - mN)) : count;
			hypergeometric_distribution<> hgd(mPlusN, aNode.supportSum, mN);
			double cutoff = pdf(hgd, count);
			double tmp_p = 0.0;

			for(uint k = min_for_k; k <= max_for_k; k++){
				double p = pdf(hgd, k);
				if(mSide == 0){
					if(p <= cutoff){
						tmp_p += p;
					}
				}else{
					tmp_p += p;
				}

			}
			mPvalueMap[tKey] = tmp_p;

		}
		tP.push_back(mPvalueMap[tKey]);
	}

	return tP;
}

void FastWY::init(const vector<vector<Event> > aTransaction, const vector<double> aY){
	mTransaction = aTransaction;
	mY = aY;
	mItemSize = mTransaction[0][0].second.size();
	mN = mTransaction.size();
	if(mTransaction[0][0].first.size() > 0) mFlagItemsetExist = 1;

	Itemset tItemset;
	vector<uint> tWild(mItemSize, MAXVAL);
	Event tEvent(tItemset, tWild);
	mWildEvent = tEvent;

	if(mTransaction.empty() || mY.empty()){
		cout << "error:Data or label is empty." << '\n';
		exit(1);
	}
	mPattern.clear();
	mTree.clear();

	mPlusN = 0;
	for(auto it : mY){
		if(it == 1) mPlusN++;
	}


	map<Event, projectDB> counter;
	map<Event, uint> dupCheck;

	if(mInterval < 0){

		for(uint i = 0; i < mN; ++i){

			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;

							Position tPosition;
							tPosition.sequence_index = i;
							tPosition.event_index = j;
							tPosition.itemset_index = k;
							counter[tEvent].push_back(tPosition);
						}
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					if(dupCheck.find(tEvent) == dupCheck.end()){
						dupCheck[tEvent] = 0;
						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						counter[tEvent].push_back(tPosition);
					}
				}
			}
			dupCheck.clear();
		}
	}else{
		for(uint i = 0; i < mN; ++i){
			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;

						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						counter[tEvent].push_back(tPosition);
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					Position tPosition;
					tPosition.sequence_index = i;
					tPosition.event_index = j;
					tPosition.itemset_index = k;
					counter[tEvent].push_back(tPosition);
				}
			}
		}
	}

	mPattern.clear();


	Node tRoot;
	mTree.push_back(tRoot);


	for(auto it = counter.begin(), end = counter.end(); it != end; ++it){る
		mPattern.push_back(it->first);
		Node tNode;
		tNode.parent = mTree.begin();
		tNode.pattern = mPattern;
		tNode.patternStr = pat2str();
		tNode.supportSum = 0;
		tNode.pdb = it->second;
		tNode.f_val = -1;
		tNode.p = 1;
		tNode.closed = true;

		uint oid = MAXVAL;

		for(uint i = 0, size = tNode.pdb.size(); i < size; ++i){
			uint id = tNode.pdb[i].sequence_index;
			if(id != oid){
				double tSup = calcSup(id, mPattern);
				tNode.supportSum += tSup;
				tNode.support.push_back(tSup);
				tNode.x.push_back(id);
			}
			oid = id;
		}

		bool tAddFlag = true;
		if(mMinsup >= 1 && mMinsup > tNode.supportSum) tAddFlag = false;
		else if(mMinsup < 1 && mMinsup > (double) tNode.supportSum / mN) tAddFlag = false;

		if(tAddFlag){
			tNode.pdb_size = sum_items_pdb(tNode.pdb);
			mTree.push_back(tNode);
			mTree.begin()->child.push_back(prev(mTree.end()));
		}
		mPattern.pop_back();
	}

}


uint FastWY::isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2){
	vector<Event> tShort_seq;
	vector<Event> tLong_seq;
	uint tSub_seq;

	if(aSeq1.size() < aSeq2.size()){
		tShort_seq = aSeq1;
		tLong_seq = aSeq2;
		tSub_seq = 1;
	}else if(aSeq1.size() > aSeq2.size()){
		tShort_seq = aSeq2;
		tLong_seq = aSeq1;
		tSub_seq = 2;
	}else if(mFlagItemsetExist){
		uint tSum1 = 0;
		uint tSum2 = 0;
		for(auto it : aSeq1){
			tSum1 += it.first.size();
		}
		for(auto it : aSeq2){
			tSum2 += it.first.size();
		}

		if(tSum1 > tSum2){
			tShort_seq = aSeq2;
			tLong_seq = aSeq1;
			tSub_seq = 2;
		}else if(tSum1 < tSum2){
			tShort_seq = aSeq1;
			tLong_seq = aSeq2;
			tSub_seq = 1;
		}else{
			return 0;
		}

	}else{
		return 0;
	}

	uint tCount = 0;
	uint diff_size = tLong_seq.size() - tShort_seq.size();

	for(uint i = 0; i <= diff_size; ++i){
		for(uint it = 0; it < tShort_seq.size(); ++it){
			if(isInclude(tLong_seq[it + i], tShort_seq[it]) || tShort_seq[it] == mWildEvent || tLong_seq[it + i] == mWildEvent){
				tCount++;
			}else{
				tCount = 0;
				break;
			}

			if(tCount == tShort_seq.size()){
				return tSub_seq;
			}
		}
	}

	return 0;
}


bool FastWY::isInclude(const Event aEvent1, const Event aEvent2){
	if(aEvent1.first.size() < aEvent2.first.size()){
		return false;
	}else if(aEvent1.second != aEvent2.second){
		return false;
	}else{
		uint i = 0;
		bool tExistFlag = false;
		for(auto itr : aEvent2.first){
			tExistFlag = false;
			if(i == aEvent1.first.size()){
				return false;
			}
			while(itr >= aEvent1.first[i]){
				if(itr == aEvent1.first[i]){
					tExistFlag = true;
					i++;
					break;
				}else{
					i++;
				}

				if(i == aEvent1.first.size()){
					return false;
				}
			}

			if(!tExistFlag){
				return false;
			}
		}
		return true;
	}

}


bool FastWY::isParent(const TreeIter aNodeIter, const TreeIter aChildIter){
	auto current = aChildIter->parent;
	while(current->pattern.size() >= aNodeIter->pattern.size()){
		if(current == aNodeIter) return true;

		current = current->parent;
	}
	return false;
}


uint FastWY::sum_items_pdb(const projectDB aPdb){
	uint tSum = 0;
	for(auto itr : aPdb){
		uint id = itr.sequence_index;
		uint j = itr.event_index;
		if(mTransaction[id][j].first.size() > 0){
			tSum += mTransaction[id][j].first.size() - 1 - itr.itemset_index;
		}else{
			tSum++;
		}
		for(++j; j < mTransaction[id].size(); ++j){
			if(mTransaction[id][j].first.size() > 0){
				tSum += mTransaction[id][j].first.size();
			}else{
				tSum++;
			}
		}
	}
	return tSum;

}


uint FastWY::sumId(const vector<uint> aIdList){
	uint tSum = 0;

	set<uint> tIdSet;

	for(auto tId: aIdList){
		tIdSet.insert(tId);
	}

	for(auto itr = tIdSet.begin();itr != tIdSet.end();++itr){
    	tSum += *itr;
	}

	return tSum;
}


uint FastWY::calcSup(uint aId, vector<Event> aPattern){
	switch(mSupMode){
		case 1:
			return calcPat(aId, aPattern);
			break;
		default:
			return 1;
			break;
	}
}


uint FastWY::calcPat(uint aId, vector<Event> aPattern){
	uint patSize = aPattern.size();
	uint k;
	uint num = 0;
	uint interval;

	if(mTransaction[aId].size() < patSize) return 0;

	for(uint i = 0; i <= mTransaction[aId].size() - patSize; i++){
		interval = 0;
		k = 0;
		for(uint j = i; j < mTransaction[aId].size(); j++){
			if(mInterval >= 0 && interval > mInterval){
				break;
			}
			if(mTransaction[aId][j] == aPattern[k] || aPattern[k] == mWildEvent){
				k++;
				interval = 0;
				if(k == patSize){
					k = 0;
					num++;
					i = j;
					break;
				}
			}else{
				interval++;
			}

		}
	}
	return num;
}


bool FastWY::calculate(Node& aNode){
	if(aNode.supportSum < mMinsup) return true;

	if(aNode.f_val == -1){
		double tLow_L;
		double tLow_U;

		if(mSide != 2){
			if(aNode.supportSum <= mPlusN){
				tLow_U = log_combination(mPlusN, aNode.supportSum) - log_combination(mN, aNode.supportSum);
			}else{
				tLow_U = 0 - log_combination(mN, mPlusN);
			}
		}

		if(mSide != 1){
			if(aNode.supportSum <= (mN - mPlusN)){
				tLow_L = log_combination(mN - mPlusN, aNode.supportSum) - log_combination(mN, aNode.supportSum);
			}else{
				tLow_L = 0 - log_combination(mN, mN - mPlusN);
			}
		}

		switch(mSide){
			case 0:
				aNode.f_val = exp(min(tLow_U, tLow_L));
				break;
			case 1:
				aNode.f_val = exp(tLow_U);
				break;
			case 2:
				aNode.f_val = exp(tLow_L);
				break;
			default:
				aNode.f_val = exp(min(tLow_U, tLow_L));
				break;
		}
	}

	vector<double> tSortL = mL;
	sort(tSortL.begin(), tSortL.end());
	if(tSortL[round(mR * mAlpha + 1)] <= aNode.f_val) return false;
	else return true;
}


bool FastWY::calculate_new(Node& aNode){
	uint oid = MAXVAL;

	for(uint i = 0, size = aNode.pdb.size(); i < size; ++i){
		uint id = aNode.pdb[i].sequence_index;
		if(oid != id){
			uint tSup = calcSup(id, aNode.pattern);
			aNode.supportSum += tSup;
			aNode.support.push_back(tSup);
			aNode.x.push_back(id);

		}
		oid = id;
	}

	double tLow_L;
	double tLow_U;

	if(mSide != 2){
		if(aNode.supportSum <= mPlusN){
			tLow_U = log_combination(mPlusN, aNode.supportSum) - log_combination(mN, aNode.supportSum);
		}else{
			tLow_U = 0 - log_combination(mN, mPlusN);
		}
	}

	if(mSide != 1){
		if(aNode.supportSum <= (mN - mPlusN)){
			tLow_L = log_combination(mN - mPlusN, aNode.supportSum) - log_combination(mN, aNode.supportSum);
		}else{
			tLow_L = 0 - log_combination(mN, mN - mPlusN);
		}
	}

	switch(mSide){
		case 0:
			aNode.f_val = exp(min(tLow_U, tLow_L));
			break;
		case 1:
			aNode.f_val = exp(tLow_U);
			break;
		case 2:
			aNode.f_val = exp(tLow_L);
			break;
		default:
			aNode.f_val = exp(min(tLow_U, tLow_L));
			break;
	}

	if(mMinsup >= 1 && mMinsup > aNode.supportSum) return false;
	else if(mMinsup < 1 && mMinsup > (double) aNode.supportSum / mN) return false;

	vector<double> tSortL = mL;
	sort(tSortL.begin(), tSortL.end());
	if(tSortL[round(mR * mAlpha + 1)] <= aNode.f_val) return false;

	return true;
}


void FastWY::project(const TreeIter aNodeIter){
	projectDB tPdb = aNodeIter->pdb;

	if(mFlagItemsetExist && mPattern.back().first.size() < mMaxItem){
		map<Event, projectDB> tCounter;
		for(uint i = 0, size = tPdb.size(); i < size; ++i){
			uint id = tPdb[i].sequence_index;
			uint j = tPdb[i].event_index;
			if(mTransaction[id][j].first.size() - 1 > tPdb[i].itemset_index){
				uint k = tPdb[i].itemset_index + 1;
				for(; k < mTransaction[id][j].first.size(); ++k){
					Event tEvent;
					vector<uint> tItemset = mPattern.back().first;
					tItemset.push_back(mTransaction[id][j].first[k]);
					tEvent.first = tItemset;
					tEvent.second = mTransaction[id][j].second;
					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;
					tPos.itemset_index = k;
					tCounter[tEvent].push_back(tPos);
				}
			}
		}

		// project: next event
		Event tSaveEvent = mPattern.back();
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.pop_back();
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
		}
		mPattern.pop_back();
		mPattern.push_back(tSaveEvent);
	}


	if(mFlagItemsetExist){
		for(uint i = 0; i < mPattern.size(); ++i){
			if(mPattern[i].first.size() < mMinItem){
				return;
			}
		}
	}

	if(mPattern.size() < mMaxpat){
		//S-Extension
		map<Event, projectDB> tCounter;
		if(mWild == 1){
			projectDB tPdbWild;

			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				Position tPosition;
				tPosition.sequence_index = tPdb[i].sequence_index;
				tPosition.event_index = tPdb[i].event_index + 1;
				tPosition.itemset_index = 0;
				if(tPosition.event_index < mTransaction[tPosition.sequence_index].size()){
					tPdbWild.push_back(tPosition);
				}
			}

			mPattern.push_back(mWildEvent);
			project_new(tPdbWild, aNodeIter);
			mPattern.pop_back();
		}

		if(mInterval < 0){
			map<Event, uint> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();

				uint j = tPdb[i].event_index + 1;

				for(; j < trsize; j++){
					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					Position tPosition;
					tPosition.sequence_index = id;
					tPosition.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							if(dupCheck.find(tEvent) == dupCheck.end()){
								dupCheck[tEvent] = 0;
								tPosition.itemset_index = k;
								tCounter[tEvent].push_back(tPosition);
							}
							tItemset.clear();
						}
					}else{
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;
							tPosition.itemset_index = k;
							tCounter[tEvent].push_back(tPosition);
						}
					}
				}
				dupCheck.clear();
			}
		}else{
			vector<Position> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();

				uint j = tPdb[i].event_index + 1;
				uint j_tmp = j;

				for(; j < trsize; j++){
					if((int) j - j_tmp > mInterval){
						break;
					}

					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							tPos.itemset_index = k;
							if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
								dupCheck.push_back(tPos);
								tCounter[tEvent].push_back(tPos);
							}
							tItemset.clear();
						}
					}else{
						tPos.itemset_index = k;
						tEvent.first = tItemset;
						if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
							dupCheck.push_back(tPos);
							tCounter[tEvent].push_back(tPos);
						}
					}
				}
			}
		}

		// project: next event
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
			mPattern.pop_back();
		}
	}

}


void FastWY::project_new(projectDB& aPdb, const TreeIter aParent){
	Node tNode;
	tNode.parent = aParent;
	tNode.pattern = mPattern;
	tNode.patternStr = pat2str();
	tNode.pdb = aPdb;
	tNode.supportSum = 0;
	tNode.f_val = -1;
	tNode.p = 1;
	tNode.closed = true;
	tNode.pdb_size = sum_items_pdb(tNode.pdb);

	TreeIter current = mTree.insert(mTree.end(), tNode);

	bool flag = calculate_new(*current);
	if(flag){

		if(mCloSpan == 1){
			if(current->parent->supportSum == current->supportSum){
				current->parent->closed = false;
			}
			mFlagCScheckEnd = 0;

			checkProjectedDB(current);
			if(mFlagCScheckEnd == 2){
				mTree.pop_back();
				return;
			}
		}

		bool tSkipFlag = false;
		if(current->pattern.size() < mMinpat){
			tSkipFlag = true;
			cout << "mMinpat skip" << endl;
		}else if(mFlagItemsetExist){
			if((current->pattern[current->pattern.size() - 1].first.size() < mMinItem)){
				tSkipFlag = true;
				cout << "mMinItem skip" << endl;
			}
		}

		if(!tSkipFlag){
			vector<double> tP = exactTest_vector(*current);

			for(uint i = 0; i < mR; ++i){
				if(tP[i] < mL[i]) mL[i] = tP[i];
			}
		}

		aParent->child.push_back(current);

		if(mFlagCScheckEnd == 0){
			project(current);
		}
	}
}

void FastWY::checkProjectedDB(const TreeIter aCurrentIter){
	uint tKey = sumId(aCurrentIter->x) + aCurrentIter->supportSum + aCurrentIter->pdb_size;

	if(mCheckPDB.find(tKey) != mCheckPDB.end()){
		for(auto itr = mCheckPDB[tKey].begin(); itr != mCheckPDB[tKey].end();){
			if(aCurrentIter->pdb_size == (*itr)->pdb_size){
				uint tWhichSub = isSubsequence(aCurrentIter->pattern, (*itr)->pattern);
				if(tWhichSub == 1){
					mFlagCScheckEnd = 2;
					return;
				}else if(tWhichSub == 2){
					aCurrentIter->child = (*itr)->child;
					(*itr)->parent->child.erase(find((*itr)->parent->child.begin(), (*itr)->parent->child.end(), *itr));
					mTree.erase(*itr);
					//
					itr = mCheckPDB[tKey].erase(itr);

					for(auto it : aCurrentIter->child){
						it->parent = aCurrentIter;
					}

					for(auto it : aCurrentIter->child){
						childPatternUpdate(it);
					}
					mFlagCScheckEnd = 1;
					return;
				}
			}

			itr++;
		}

	}

	mCheckPDB[tKey].push_back(aCurrentIter);
	return;
}

void FastWY::childPatternUpdate(const TreeIter aChildIter){
	Event tChildLast = aChildIter->pattern.back();

	Event tParentLast = mPattern.back();
	bool tIE = false;
	if(tChildLast.first.size() > tParentLast.first.size()){
		mPattern.pop_back();
		tIE = true;
	}
	mPattern.push_back(tChildLast);
	aChildIter->pattern = mPattern;
	aChildIter->patternStr = pat2str();
	for(auto it : aChildIter->child){
		childPatternUpdate(it);
	}
	mPattern.pop_back();
	if(tIE) mPattern.push_back(tParentLast);
}

void FastWY::search_tree(const TreeIter aNodeIter){

	bool flag = calculate(*aNodeIter);
	if(flag){

		bool tSkipFlag = false;
		if(aNodeIter->pattern.size() < mMinpat){
			tSkipFlag = true;
			cout << "mMinpat skip" << endl;
		}else if(mFlagItemsetExist){
			if((aNodeIter->pattern[aNodeIter->pattern.size() - 1].first.size() < mMinItem)){
				tSkipFlag = true;
				cout << "mMinItem skip" << endl;
			}
		}

		if(!tSkipFlag){
			vector<double> tP = exactTest_vector(*aNodeIter);

			for(uint i = 0; i < mR; ++i){
				if(tP[i] < mL[i]) mL[i] = tP[i];
			}
		}

		project(aNodeIter);

	}
}

void FastWY::search_tree_final(const TreeIter aNodeIter){
	if(aNodeIter->f_val <= mDelta){
		bool tSkipFlag = false;
		if(aNodeIter->pattern.size() < mMinpat){
			tSkipFlag = true;
			cout << "mMinpat skip" << endl;
		}else if(mFlagItemsetExist){
			if((aNodeIter->pattern[aNodeIter->pattern.size() - 1].first.size() < mMinItem)){
				tSkipFlag = true;
				cout << "mMinItem skip" << endl;
			}
		}
		if(!aNodeIter->closed){
			tSkipFlag = true;
		}

		if(!tSkipFlag){
			aNodeIter->p = exactTest(*aNodeIter);

			if(aNodeIter->p <= mDelta){
				uint tKey = sumId(aNodeIter->x) + aNodeIter->supportSum;
				bool tAddFlag = true;
				if(mCloSpan == 1 && mSignificant.find(tKey) != mSignificant.end()){
					for(auto itr = mSignificant[tKey].begin(); itr != mSignificant[tKey].end();){
						if(aNodeIter->supportSum == (*itr)->supportSum){
							uint tWhichSub = isSubsequence(aNodeIter->pattern, (*itr)->pattern);
							if(tWhichSub == 1){
								tAddFlag = false;
								break;
							}else if(tWhichSub == 2){
								itr = mSignificant[tKey].erase(itr);
								continue;
							}
						}
						itr++;
					}
				}

				if(tAddFlag){
					mSignificant[tKey].push_back(aNodeIter);
				}

			}
		}

		if(aNodeIter->child.size() != 0){
			for(auto it : aNodeIter->child){
				search_tree_final(it);
			}
		}
	}
}



void FastWY::main(){

	mt19937 rand_src(1);

	TreeIter root = mTree.begin();

	for(uint rep = 0; rep < mR; ++rep){

		mTempY.push_back(mY);
		for(uint j = 0; j < mN; ++j){
			uint k = j + (rand_src() % (mN - j));
			swap(mTempY[rep][j], mTempY[rep][k]);
		}

		mL.push_back(mAlpha);
	}

	for(auto it : root->child){
		mPattern = (it->pattern);
		search_tree(it);
	}

	sort(mL.begin(), mL.end());

	double tMax = mL[round(mR * mAlpha + 1)];

	int it = round(mR * mAlpha + 1) - 1;

	while(it >= 0){
		if(mL[it] < tMax){
			mDelta = mL[it];
			break;
		}else it--;
	}

	cout << "delta:" << mDelta << endl;
	cout << "FWER = " << (double)it / mR << endl;

	for(auto it : root->child){
		search_tree_final(it);
	}

}

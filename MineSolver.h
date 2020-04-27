#ifndef MINE_MINESOLVER_H
#define MINE_MINESOLVER_H

#include <set>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

struct MineSolver {
	int SUMLIMIT = 10000000;
	int DEEPDFSLIMIT = 1000;

	int n, m, mine;
	vector <vector <long double> > mineProb;

	MineSolver (int n, int m, int mine) : n(n), m(m), mine(mine) {
		mineProb = std::vector <vector <long double> > (n * m + 1, vector <long double> (mine + 2, 0));
		mineProb[0][0] = 1.0;
		for (int i = 0; i < n * m; i++) {
			for (int j = 0; j <= mine; j++) {
				mineProb[i + 1][j] += mineProb[i][j];
				mineProb[i + 1][j + 1] += mineProb[i][j];
			}
		}
	}

	void SimpleDetect(vector <vector <int> > &board) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (board[i][j] > 0) {
					int s = 0;
					for (int ii = max(0, i - 1); ii < min(n, i + 2); ii++) {
						for (int jj = max(0, j - 1); jj < min(m, j + 2); jj++) {
							s += (board[ii][jj] < 0);
						}
					}
					if (s == board[i][j]) {
						for (int ii = max(0, i - 1); ii < min(n, i + 2); ii++) {
							for (int jj = max(0, j - 1); jj < min(m, j + 2); jj++) {
								board[ii][jj] -= (board[ii][jj] == -1);
							}
						}
					}
				}
			}
		}
	}

	bool flag;
	void DFS(vector <pair <int, pair <int, int> > > &cnt, vector <vector <int> > &edges, int n, int i,
	         vector <int> &s, vector <int> &cur, int &sum, int mineMin, int mineMax,
	         std::map <int, long double> &mineCnt, vector <std::map <int, long double> > &perBoardCnt,
	         vector <vector <int> > &allStatus, bool saveAllStatus) {
		if (sum == SUMLIMIT) {
			flag = false;
			return;
		}
		if (mineMax < 0 || mineMin > (n - i)) {
			return;
		}
		for (int j = 0; j < cnt.size(); j++) {
			if (cnt[j].first < 0 || cnt[j].second.first + cnt[j].second.second < cnt[j].first) {
				return;
			}
		}
		if (i == n) {
			sum++;
			int mineCur = 0;
			for (int i = 0; i < s.size(); i++) {
				s[i] += cur[i];
				mineCur += cur[i];
			}
			mineCnt[mineCur]++;
			for (int i = 0; i < s.size(); i++) {
				if (cur[i]) {
					perBoardCnt[i][mineCur]++;
				}
			}
			if (saveAllStatus) {
				allStatus.push_back(cur);
				for (int i = 0; i < allStatus.back().size(); i++) {
					if (allStatus.back()[i] == 1) {
						allStatus.back()[i] = -2;
					}
				}
			}
			return;
		}
		for (int j = 0; j < edges[i].size(); j++) {
			cnt[edges[i][j]].second.first--;
		}
		cur[i] = 0;
		DFS(cnt, edges, n, i + 1, s, cur, sum, mineMin, mineMax, mineCnt, perBoardCnt, allStatus, saveAllStatus);
		for (int j = 0; j < edges[i].size(); j++) {
			cnt[edges[i][j]].second.first++;
		}
		cur[i] = 1;
		for (int j = 0; j < edges[i].size(); j++) {
			cnt[edges[i][j]].first--;
			cnt[edges[i][j]].second.first--;
		}
		DFS(cnt, edges, n, i + 1, s, cur, sum, mineMin - 1, mineMax - 1, mineCnt, perBoardCnt, allStatus, saveAllStatus);
		cur[i] = 0;
		for (int j = 0; j < edges[i].size(); j++) {
			cnt[edges[i][j]].first++;
			cnt[edges[i][j]].second.first++;
		}
	}

	void DFSDetect(vector <vector <int> > &board, vector <pair <int, int> > &points, vector <pair <int, int> > &safeCells, vector <pair <int, int> > &unknownCells, vector <vector <double> > &prob,
	               int mineMin, int mineMax, std::map <int, long double> &mineCnt, vector <std::map <int, long double> > &perBoardCnt, vector <vector <int> > &allStatus, bool saveAllStatus) {
		set <pair <int, int> > unknownSet;
		for (int i = 0; i < (int)unknownCells.size(); i++) {
			unknownSet.insert(unknownCells[i]);
		}
		vector <pair <int, pair <int, int> > > cnt;
		vector <vector <int> > edges;
		for (int i = 0; i < (int)safeCells.size(); i++) {
			pair <int, int> &it = safeCells[i];
			cnt.push_back(make_pair(board[it.first][it.second], make_pair(0, 0)));
			for (int ii = max(0, it.first - 1); ii < min(n, it.first + 2); ii++) {
				for (int jj = max(0, it.second - 1); jj < min(m, it.second + 2); jj++) {
					cnt.back().first -= (board[ii][jj] == -2);
					if (board[ii][jj] == -1) {
						if (unknownSet.find(make_pair(ii, jj)) == unknownSet.end()) {
							cnt.back().second.second++;
						} else {
							cnt.back().second.first++;
						}
					}
				}
			}
		}
		for (int i = 0; i < (int)unknownCells.size(); i++) {
			edges.push_back(vector<int>());
			for (int j = 0; j < (int)safeCells.size(); j++) {
				if (abs(safeCells[j].first - unknownCells[i].first) <= 1 && abs(safeCells[j].second - unknownCells[i].second) <= 1) {
					edges[i].push_back(j);
				}
			}
		}
		int cells = unknownCells.size();
		int sum = 0;
		vector <int> s(cells, 0), cur(cells, 0);
		mineCnt.clear();
		perBoardCnt.clear();
		perBoardCnt.resize(cells);
		DFS(cnt, edges, cells, 0, s, cur, sum, mineMin, mineMax, mineCnt, perBoardCnt, allStatus, saveAllStatus);
		for (int i = 0; i < cells; i++) {
			if (s[i] == 0) {
				points.push_back(unknownCells[i]);
			}
			else if (s[i] == sum) {
				board[unknownCells[i].first][unknownCells[i].second] = -2;
			} else {
				prob[unknownCells[i].first][unknownCells[i].second] = (double)s[i] / sum;
			}
		}
		if (saveAllStatus) {
			for (int i = 0; i < allStatus.size(); i++) {
				vector <int> &curStatus = allStatus[i];
				for (int j = 0; j < curStatus.size(); j++) {
					if (curStatus[j] == 0) {
						for (int ii = max(0, unknownCells[j].first - 1); ii < min(n, unknownCells[j].first + 2); ii++) {
							for (int jj = max(0, unknownCells[j].second - 1); jj < min(m, unknownCells[j].second + 2); jj++) {
								curStatus[j] += (board[ii][jj] == -2);
							}
						}
						for (int k = 0; k < curStatus.size(); k++) {
							curStatus[j] += (curStatus[k] == -2 &&
							                 abs(unknownCells[j].first - unknownCells[k].first) <= 1 &&
							                 abs(unknownCells[j].second - unknownCells[k].second) <= 1);
						}
					}
				}
			}
		}
	}

	void Floodfill(vector <vector <int> > &board, set <pair <int, int> > &use, vector <pair <int, int> > &safeCells, vector <pair <int, int> > &unkownCells, int i, int j) {
		use.insert(make_pair(i, j));
		if (board[i][j] == -1) {
			unkownCells.push_back(make_pair(i, j));
		} else {
			safeCells.push_back(make_pair(i, j));
		}
		for (int ii = max(0, i - 1); ii < min(n, i + 2); ii++) {
			for (int jj = max(0, j - 1); jj < min(m, j + 2); jj++) {
				if (board[ii][jj] != -2 && use.find(make_pair(ii, jj)) == use.end() &&
				    (board[ii][jj] * board[i][j] < 0)) {
					Floodfill(board, use, safeCells, unkownCells, ii, jj);
				}
			}
		}
	}

	void StatusSplit(vector <vector <int> > &allStatus, vector <int> ids, vector <vector <int> > &idsList) {
		if (ids.size() == 0) {
			return;
		}
		for (int i = 0; i < allStatus[0].size(); i++) {
			bool allSame = true, noMine = true;
			for (int id : ids) {
				allSame &= (allStatus[id][i] == allStatus[ids[0]][i]);
				noMine &= (allStatus[id][i] >= 0);
			}
			if (noMine && !allSame) {
				std::map <int, vector <int> > idsDict;
				for (int id : ids) {
					idsDict[allStatus[id][i]].push_back(id);
				}
				for (auto spids : idsDict) {
					StatusSplit(allStatus, spids.second, idsList);
				}
				return;
			}
		}
		idsList.push_back(ids);
	}

	double DeepDFS(vector <vector <int> > &allStatus, vector <int> &ids, vector <double> &prob, std::map<vector <int>, double> &f) {
		if (f.find(ids) != f.end()) {
			return f[ids];
		}
		if (ids.size() == 1) {
			for (int sel = 0; sel < prob.size(); sel++) {
				prob[sel] = (allStatus[ids[0]][sel] < 0);
			}
			return f[ids] = 1.0;
		}
		double ret = 0.0;
		vector <pair <int, int> > order;
		for (int sel = 0; sel < prob.size(); sel++) {
			int mine = 0;
			for (int i : ids) {
				mine += (allStatus[i][sel] < 0);
			}
			order.push_back(make_pair(mine, sel));
		}
		for (auto it : order) {
			int sel = it.second;
			vector <vector <int> > idsList;
			vector <int> cur;
			for (int i : ids) {
				if (allStatus[i][sel] >= 0) {
					cur.push_back(i);
				}
			}
			if (cur.size() == 0) {
				prob[sel] = 0;
				continue;
			}
			if (cur.size() == ids.size()) {
				prob[sel] = 1;
				continue;
			}
			if ((double)cur.size() / ids.size() < ret - 1e-8) {
				prob[sel] = 0;
				continue;
			}
			StatusSplit(allStatus, cur, idsList);
			vector <double> p = vector<double>(prob.size(), 0);
			prob[sel] = 0.0;
			for (int i = 0; i < idsList.size(); i++) {
				prob[sel] += DeepDFS(allStatus, idsList[i], p, f) * idsList[i].size() / ids.size();
			}
			ret = max(ret, prob[sel]);
		}
		return f[ids] = ret;
	}

	bool DFSDetect(vector <vector <int> > &board, vector <pair <int, int> > &points, vector <vector <double> > &prob) {
		int unkownS = 0, surMine = mine;
		int rangeMin, rangeMax;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				unkownS += (board[i][j] == -1);
				surMine -= (board[i][j] == -2);
			}
		}
		std::map <int, long double> mineCnt;
		vector <std::map <int, long double> > perBoardCnt;

		vector <vector <pair <int, int> > > safeCellsList, unkownCellsList;
		vector <std::map <int, long double> > mineCntList;
		vector <vector <std::map <int, long double> > > perBoardCntList;

		vector <vector <int> > ori = board;
		vector <vector <int> > allStatus;

		set <pair <int, int> > use;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (board[i][j] != -2 && use.find(make_pair(i, j)) == use.end()) {
					vector <pair <int, int> > safeCells, unkownCells;
					Floodfill(board, use, safeCells, unkownCells, i, j);
					if (safeCells.size() >= 1 && unkownCells.size() >= 1) {
						DFSDetect(board, points, safeCells, unkownCells, prob, 0, surMine, mineCnt, perBoardCnt, allStatus, false);
						if (points.size() > 0) {
							return true;
						}
						safeCellsList.push_back(safeCells);
						unkownCellsList.push_back(unkownCells);
						mineCntList.push_back(mineCnt);
						perBoardCntList.push_back(perBoardCnt);
					}
				}
			}
		}

		int last = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				last += (board[i][j] == -1 && prob[i][j] < -0.5);
			}
		}

		bool doDeepDFS = false;
		for (int i = 0; i < safeCellsList.size(); i++) {
			vector <long double> f = mineProb[last];
			for (int j = 0; j < safeCellsList.size(); j++) {
				if (i != j) {
					for (int k = surMine; k >= 0; k--) {
						long double s = 0.0;
						for (auto it : mineCntList[j]) {
							if (k >= it.first) {
								s += f[k - it.first] * it.second;
							}
						}
						f[k] = s;
					}
				}
			}
			long double s = 0.0;
			for (auto it : mineCntList[i]) {
				if (it.first <= surMine) {
					s += it.second * f[surMine - it.first];
				}
			}

			if (s <= DEEPDFSLIMIT) {
				doDeepDFS = true;
				break;
			}

			for (int j = 0; j < unkownCellsList[i].size(); j++) {
				long double curS = 0.0;
				for (auto it : perBoardCntList[i][j]) {
					if (it.first <= surMine) {
						curS += it.second * f[surMine - it.first];
					}
				}

				prob[unkownCellsList[i][j].first][unkownCellsList[i][j].second] = curS / s;
			}
		}

		if (doDeepDFS) {
			vector<pair<int, int> > safeCells, unkownCells;
			surMine = mine;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					if (board[i][j] == -1) {
						unkownCells.push_back(make_pair(i, j));
					} else if (board[i][j] > 0) {
						safeCells.push_back(make_pair(i, j));
					}
					surMine -= (board[i][j] == -2);
				}
			}

			DFSDetect(board, points, safeCells, unkownCells, prob, surMine, surMine, mineCnt, perBoardCnt, allStatus, true);
			for (int i = 0; i < unkownCells.size(); i++) {
				bool noMine = true;
				for (int j = 0; j < allStatus.size(); j++) {
					if (allStatus[j][i] < 0) {
						noMine = false;
						break;
					}
				}
				if (noMine) {
					points.push_back(unkownCells[i]);
				}
			}
			if (points.size() > 0) {
				return true;
			}

			vector <double> deepProb = vector <double> (unkownCells.size());
			vector <int> ids;
			for (int i = 0; i < allStatus.size(); i++) {
				ids.push_back(i);
			}
			std::map <vector <int>, double> f;
			DeepDFS(allStatus, ids, deepProb, f);
			for (int i = 0; i < unkownCells.size(); i++) {
				prob[unkownCells[i].first][unkownCells[i].second] = 1 - deepProb[i];
			}
		} else {
			long double psum = 0.0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					if (board[i][j] == -1 && prob[i][j] > -0.5) {
						psum += prob[i][j];
					}
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					if (board[i][j] == -1 && prob[i][j] < -0.5) {
						prob[i][j] = (surMine - psum) / last;
					}
				}
			}
		}

		return false;
	}

	bool GetSafePoints(vector <vector <int> > &board, vector <pair <int, int> > &points) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (board[i][j] >= 0) {
					int s = 0;
					for (int ii = max(0, i - 1); ii < min(n, i + 2); ii++) {
						for (int jj = max(0, j - 1); jj < min(m, j + 2); jj++) {
							s += (board[ii][jj] == -2);
						}
					}
					if (s == board[i][j]) {
						for (int ii = max(0, i - 1); ii < min(n, i + 2); ii++) {
							for (int jj = max(0, j - 1); jj < min(m, j + 2); jj++) {
								if (board[ii][jj] == -1) {
									points.push_back(make_pair(ii, jj));
								}
							}
						}
					}
				}
			}
		}
		return points.size() > 0;
	}

	int GetValue(int i, int j, vector <vector <int> > &board) {
		if ((i == 0 || i == n - 1) && (j == 0 || j == m - 1)) {
			return -1000;
		} else {
			int s = 0;
			for (int ii = max(0, i - 2); ii < min(n, i + 3); ii++) {
				for (int jj = max(0, j - 2); jj < min(m, j + 3); jj++) {
					s += (board[ii][jj] >= 0);
				}
			}
			return -s;
		}
	}

	bool GetCLKPoints(vector <vector <int> > board, vector <pair <int, int> > &points) {
		flag = true;
		points.clear();
		SimpleDetect(board);
		if (GetSafePoints(board, points)) {
			return true;
		}
		vector <vector <double> > p(n, vector <double>(m, -1));
		if (DFSDetect(board, points, p)) {
			return flag;
		}
		if (board[0][0] == -1) {
			points.push_back(make_pair(0, 0));
			return true;
		}
		double minP = 1e100, totalNone = 0, totalMine = 0;
		int si = -1, sj = -1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (board[i][j] == -1 &&
				    (p[i][j] < minP - 1e-8 ||
				     (fabs(p[i][j] - minP) < 1e-8 && GetValue(i, j, board) < GetValue(si, sj, board)))) {
					minP = p[i][j];
					si = i;
					sj = j;
				}
			}
		}
		points.push_back(make_pair(si, sj));
		return (minP < 1e-8) && flag;
	}
};
#endif //MINE_MINESOLVER_H

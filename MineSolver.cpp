#include <ctime>
#include "MineSolver.h"

void GenMap(int n, int m, int mine, vector <vector <int> > &map, vector <vector <int> > &board) {
	vector <pair <int, int> > points;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (i != 0 || j != 0) {
				points.push_back(make_pair(i, j));
			}
		}
	}

	map.clear();
	map = std::vector <vector <int> >(n, vector <int> (m, 0));
	board.clear();
	board = std::vector <vector <int> >(n, vector <int> (m, -1));
	random_shuffle(points.begin(), points.end());
	points.erase(points.begin() + mine, points.end());
	for (auto point : points) {
		map[point.first][point.second] = -2;
	}
	for (auto point : points) {
		for (int i = max(0, point.first - 1); i < min(n, point.first + 2); i++) {
			for (int j = max(0, point.second - 1); j < min(m, point.second + 2); j++) {
				map[i][j] += (map[i][j] >= 0);
			}
		}
	}
}

int main() {
	srand(time(NULL));
    int n = 16, m = 30, mine = 99, right = 0;
    MineSolver solver(n, m, mine);
    vector <vector <int> > board, map;
    for (int times = 0; times < 100000; times++) {
        GenMap(n, m, mine, map, board);
        set <pair <int, int> > opens;
        bool succ = true;
        while (opens.size() < n * m - mine && succ) {
            vector <pair <int, int> > points;
            solver.GetCLKPoints(board, points);
            for (auto point : points) {
                opens.insert(point);
                board[point.first][point.second] = map[point.first][point.second];
                succ &= (map[point.first][point.second] != -2);
            }
        }
        right += succ;

        if (times % 100 == 99) {
            printf("%d / %d, %f\r", right, (times + 1), (double)right / (double)times);
            fflush(stdout);
        }
    }
    printf("right = %d\n", right);
    return 0;
}

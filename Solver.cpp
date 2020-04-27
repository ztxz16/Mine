#include "MineSolver.h"

int main() {
	freopen("board.txt", "r", stdin);
	freopen("out.txt", "w", stdout);
	int n, m, mine;
	scanf("%d %d %d", &n, &m, &mine);
	vector <vector <int> > board(n, vector <int> (m, -1));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			scanf("%d", &board[i][j]);
		}
	}
	vector <pair <int, int> > points;
	MineSolver solver(n, m, mine);
	solver.GetCLKPoints(board, points);
	for (auto it : points) {
		printf("%d %d\n", it.first, it.second);
	}
	fclose(stdin);
	fclose(stdout);
	return 0;
}
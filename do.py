import numpy as np;
from PIL import ImageGrab;
import win32api as api;
import win32gui as gui;
import win32con as con;
import time;
import os;

cell = 16;

def getBoard(x0, y0, x1, y1):
	img = np.array(ImageGrab.grab((x0, y0, x1, y1)));
	m, n = (x1 - x0) // 16, (y1 - y0) // 16;
	mine = 99;
	if (m == 9 and n == 9):
		mine = 10;
	if (m == 16 and n == 16):
		mine = 40;
	board = [[-1 for i in range(m)] for j in range(n)];
	map = {147 : -1, 527 : 1, 215 : 0, 722 : 2, 125 : 3, 452 : 4, 730 : 5, 668 : -2, 749 : -2, 377 : 6, 758 : 7, 346 : 8, 990 : 100};
	
	succ = False;
	for i in range(n):
		for j in range(m):
			cur = img[i * cell : (i + 1) * cell, j * cell : (j + 1) * cell, 0 : 1];
			s = np.sum(cur) % 999;
			if s not in map:
				print(i, j, s);
				exit(0);
			else:
				board[i][j] = map[s];
			if (board[i][j] == -2):
				return [];
			if (board[i][j] == 100):
				succ = True;
	if (succ):
		return [(-1, -1)];
	
	with open("board.txt", "w") as fo:
		fo.write(str(n) + " " + str(m) + " " + str(mine) + "\n");
		for i in range(n):
			for j in range(m):
				fo.write(str(board[i][j]) + " ");
			fo.write("\n");
	os.system("Solver.exe");
	ret = [];
	with open("out.txt", "r") as fi:
		for line in fi.readlines():
			cur = line.split(' ');
			ret.append((int(cur[0]), int(cur[1])));
	return ret;
				
handle = gui.FindWindow(None, "扫雷")
gui.SetForegroundWindow(handle);

x0, y0, x1, y1 = gui.GetWindowRect(handle);
x0, y0 = x0 + 15, y0 + 101;
#x1, y1 = x1 - 10, y1 - 40;

succ = 0;
fail = 0;

while True:
	ret = getBoard(x0, y0, x1, y1);
	if (len(ret) == 0):
		fail += 1;
		api.keybd_event(113, 0, 0, 0);
		api.keybd_event(113, 0, con.KEYEVENTF_KEYUP, 0);
		time.sleep(0.5);
		print(succ, " / ", succ + fail, "(", succ / (succ + fail), ")");
	elif (ret[0][0] == -1):
		succ += 1;
		api.keybd_event(113, 0, 0, 0);
		api.keybd_event(113, 0, con.KEYEVENTF_KEYUP, 0);
		time.sleep(0.5);
		print(succ, " / ", succ + fail, "(", succ / (succ + fail), ")");
	for i, j in ret:
		api.SetCursorPos((x0 + j * cell + 5, y0 + i * cell + 5));
		api.mouse_event(con.MOUSEEVENTF_LEFTDOWN, 0, 0, 0, 0);
		api.mouse_event(con.MOUSEEVENTF_LEFTUP, 0, 0, 0, 0);
		api.SetCursorPos((0, 0));
		#time.sleep(0.02);
	
MineSolver.h是核心文件，里面包装了一个MineSolver类，用于传入地图，返回一系列点

MineSolver.cpp是模拟测试程序

Solver.cpp是一个简单的demo，用于文件读取局面、输出落点

do.py是python前端，打开winmine.exe之后运行do.py就可以自动扫雷

do.py的依赖：

pip install pillow

pip install pywin32

xp规则下，高级模拟胜率为40.07%

win7规则下，高级模拟胜率为52.98%
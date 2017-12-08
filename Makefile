all: find-expr

find-expr: find-expr.cpp
	g++ -O3 -lpthread -o find-expr find-expr.cpp

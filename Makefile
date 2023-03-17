all: cgal-instr valid

cgal-instr: main.cpp
	afl-clang-lto++ -lgmp -lmpfr -lgmpxx -o cgal-instr -static main.cpp

valid: valid.cpp
	clang++ -o valid valid.cpp


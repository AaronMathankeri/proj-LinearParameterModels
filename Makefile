CC=g++
FLAGS=-std=c++11 #-mkl -I ./include/ -Wl,-rpath,${MKLROOT}/lib 
src = $(wildcard ./src/*.cpp)

all: compile run

compile:
	${CC} ${FLAGS} ${src} -o bin/a.out
run:
	./bin/a.out

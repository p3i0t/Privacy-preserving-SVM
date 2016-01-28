all: compile

compile: Project.cpp
	icc  -pg -mkl -static-intel -o project Project.cpp

clean:
	$(RM) a.out



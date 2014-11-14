#IDIR =/usr/include/
code_repo_path=/home/t30/ger/ga79moz/grad_research_phd/project/code_repository/

CC=g++
CFLAGS  =  -ggdb -I/sw/include -I${code_repo_path}
LDFLAGS =  -L/sw/lib -lgsl -lm -lgslcblas



DEPS = void_numerics.h void.h ${code_repo_path}bren_lib.h
OBJ  = void_numerics_driver.obj void_numerics.obj 

%.obj: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

void_numerics.x: $(OBJ)
	g++ -o $@ $^ $(LDFLAGS)

#void_numerics.db: $(OBJ)
#	g++ -g -o $@ $^ $(LDFLAGS)

#---- to compile with gprofiler data:
#---- g++ -pg -ggdb void_numerics_driver.cpp void_numerics.cpp -o void_numerics.x   -L/sw/lib -lgsl -lm -lgslcblas  -I/opt/local/include  -I/home/t30/ger/ga79moz/grad_research_phd/project/code_repository/


ALL_OBJS = void_numerics_driver.obj  void_numerics.obj
clean:
	$(RM) $(ALL_OBJS)

debug:
	g++ -g -ggdb void_numerics_driver.cpp void_numerics.cpp -o void_numerics.x $(LDFLAGS) -I/opt/local/include -I${code_repo_path}




#IDIR =/usr/include/
CC=g++
CFLAGS  =  -ggdb -I/sw/include -I/home/b/Brendan.Osberg/grad_research_phd/project/code_repository/
LDFLAGS =  -L/sw/lib -lgsl -lm -lgslcblas



DEPS = void_numerics.h void.h /home/b/Brendan.Osberg/grad_research_phd/project/code_repository/bren_lib.h
OBJ  = void_numerics_driver.obj void_numerics.obj 

%.obj: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

void_numerics.x: $(OBJ)
	g++ -o $@ $^ $(LDFLAGS)

#void_numerics.db: $(OBJ)
#	g++ -g -o $@ $^ $(LDFLAGS)

ALL_OBJS = void_numerics_driver.obj  void_numerics.obj
clean:
	$(RM) $(ALL_OBJS)

debug:
	g++ -ggdb void_numerics_driver.cpp void_numerics.cpp $(LDFLAGS) -I/opt/local/include




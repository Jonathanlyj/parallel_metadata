CC = mpicc
CFLAGS = -O0 -g 
INCLUDES = -I$(PNETCDF_DIR)/include
LFLAGS = -L$(PNETCDF_DIR)/lib
LIBS = -lpnetcdf

SRCS = baseline_ex1.c baseline_ncx.c baseline_ex2.c
OBJS = $(SRCS:.c=.o)



all: baseline_ex1 baseline_ex2

baseline_ex1: baseline_ex1.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

baseline_ex2: baseline_ex2.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS)baseline_ex1 baseline_ex2


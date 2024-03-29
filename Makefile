CC = mpicc
PNETCDF_DIR =  /global/homes/y/yll6162/pnetcdf-source/pnetcdf-install
CFLAGS = -O0 -g
INCLUDES = -I$(PNETCDF_DIR)/include
LFLAGS = -L$(PNETCDF_DIR)/lib
LIBS = -lpnetcdf

SRCS = baseline_ex1.c baseline_ncx.c app_baseline_test_all.c baseline_test.c lib_level_baseline_test_shared.c lib_level_baseline_test_read.c lib_level_baseline_test_dup_name.c lib_baseline_test_all.c
OBJS = $(SRCS:.c=.o)



all: baseline_test baseline_ex1 app_baseline_test_all pnc_consist_check lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all

baseline_test: baseline_test.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)
baseline_ex1: baseline_ex1.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

app_baseline_test_all: app_baseline_test_all.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

pnc_consist_check: pnc_consist_check.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_shared: lib_level_baseline_test_shared.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_read: lib_level_baseline_test_read.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_dup_name: lib_level_baseline_test_dup_name.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_baseline_test_all: lib_baseline_test_all.o baseline_ncx.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS)baseline_ex1 app_baseline_test_all pnc_consist_check lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all


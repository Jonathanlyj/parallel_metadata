# PNETCDF_DIR=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-install
# PNETCDF_DIR=/files2/scratch/yll6162/pnetcdf/PnetCDF-install
# PNETCDF_DIR=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-sort-install
# PNETCDF_DIR=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-format-install
PNETCDF_DIR=/homes/yll6162/PnetCDF_meta/PnetCDF-install

CC = mpicc
H5CC = h5pcc
CFLAGS = -O0 -g -I$(PNETCDF_DIR)/include 
INCLUDES = -I$(PNETCDF_DIR)/include
LFLAGS = -L$(PNETCDF_DIR)/lib
LIBS = -lpnetcdf

SRCS = baseline_ex1.c baseline_ncx_app.c baseline_ncx_lib.c app_baseline_test_all.c baseline_test.c lib_level_baseline_test_shared.c lib_level_baseline_test_read.c lib_level_baseline_test_dup_name.c lib_baseline_test_all.c h5_baseline_test_all.c
OBJS = $(SRCS:.c=.o)

LIB_PROGRAMS = lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all

all: baseline_test baseline_ex1 app_baseline_test_all pnc_consist_check $(LIB_PROGRAMS) h5_baseline_test_all

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

baseline_test: baseline_test.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)
baseline_ex1: baseline_ex1.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

app_baseline_test_all: app_baseline_test_all.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

pnc_consist_check: pnc_consist_check.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_shared: lib_level_baseline_test_shared.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_read: lib_level_baseline_test_read.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_dup_name: lib_level_baseline_test_dup_name.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_baseline_test_all: lib_baseline_test_all.o baseline_ncx_lib.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

new_format_create_simple: new_format_create_simple.o $(PNETCDF_DIR)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

new_format_create_diff: new_format_create_diff.o $(PNETCDF_DIR)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

new_format_open: new_format_open.o $(PNETCDF_DIR)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

h5_baseline_test_all: h5_baseline_test_all.o
	$(H5CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

h5_%.o: h5_%.c
	$(H5CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

copy:
	cp ./lib_level_baseline_old.nc ./lib_level_baseline.nc
clean:
	rm -f $(OBJS)baseline_ex1 app_baseline_test_all pnc_consist_check lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all h5_baseline_test_all new_format_create_simple
base: baseline_test baseline_ex1 app_baseline_test_all pnc_consist_check
lib: $(LIB_PROGRAMS)
h5: h5_baseline_test_all
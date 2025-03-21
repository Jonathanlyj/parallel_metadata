PNETCDF_DIR=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-install
# PNETCDF_DIR=/files2/scratch/yll6162/pnetcdf/PnetCDF-install
# PNETCDF_DIR=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-sort-install
PNETCDF_DIR_FORMAT=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-format-install
# PNETCDF_DIR=/homes/yll6162/PnetCDF_meta/PnetCDF-install
PNETCDF_DIR_LIBBASE=/files2/scratch/yll6162/PnetCDF-meta/PnetCDF-lib-install

CC = mpicc
H5CC = h5pcc

CFLAGS = -O2
INCLUDES = -I$(PNETCDF_DIR)/include
INCLUDES_LIBBASE = -I$(PNETCDF_DIR_LIBBASE)/include
INCLUDES_FORMAT = -I$(PNETCDF_DIR_FORMAT)/include

LFLAGS = -L$(PNETCDF_DIR)/lib
LFLAGS_LIBBASE = -L$(PNETCDF_DIR_LIBBASE)/lib
LFLAGS_FORMAT = -L$(PNETCDF_DIR_FORMAT)/lib
LIBS = -lpnetcdf

SRCS = baseline_ex1.c baseline_ncx_app.c baseline_ncx_lib.c app_baseline_test_all.c baseline_test.c lib_level_baseline_test_shared.c lib_level_baseline_test_read.c lib_level_baseline_test_dup_name.c lib_baseline_test_all.c h5_baseline_test_all.c
OBJS = $(SRCS:.c=.o)

LIB_PROGRAMS = lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all

all: baseline_test baseline_ex1 app_baseline_test_all pnc_consist_check $(LIB_PROGRAMS) h5_baseline_test_all

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

lib_%.o: lib_%.c
	$(CC) $(CFLAGS) $(INCLUDES_LIBBASE) -c $< -o $@

new_format_%.o: new_format_%.c
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) -c $< -o $@



baseline_test: baseline_test.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)
baseline_ex1: baseline_ex1.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

app_baseline_test_all: app_baseline_test_all.o baseline_ncx_app.o $(PNETCDF_DIR)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

app_baseline_read_test_all: app_baseline_read_test_all.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)
save_input_test_all: save_input_test_all.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

pnc_consist_check: pnc_consist_check.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_level_baseline_test_shared: lib_level_baseline_test_shared.o $(PNETCDF_DIR_LIBBASE)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_LIBBASE) $(LFLAGS_LIBBASE) -o $@ $^ $(LIBS)

lib_level_baseline_test_read: lib_level_baseline_test_read.o $(PNETCDF_DIR_LIBBASE)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_LIBBASE) $(LFLAGS_LIBBASE) -o $@ $^ $(LIBS)

lib_level_baseline_test_dup_name: lib_level_baseline_test_dup_name.o $(PNETCDF_DIR_LIBBASE)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_LIBBASE) $(LFLAGS_LIBBASE) -o $@ $^ $(LIBS)

lib_baseline_test_all: lib_baseline_test_all.o baseline_ncx_lib.o $(PNETCDF_DIR_LIBBASE)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_LIBBASE) $(LFLAGS_LIBBASE) -o $@ $^ $(LIBS)

new_format_create_simple: new_format_create_simple.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) -o $@ $^ $(LIBS)

new_format_create_large: new_format_create_large.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) -o $@ $^ $(LIBS)

new_format_create_conflict_blk: new_format_create_conflict_blk.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) -o $@ $^ $(LIBS)

new_format_create_diff: new_format_create_diff.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) -o $@ $^ $(LIBS)

new_format_open: new_format_open.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) -o $@ $^ $(LIBS)

new_format_test_all: new_format_test_all.o baseline_ncx_app.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) $(LIBS) -o $@ $^ 

new_format_read_test_all: new_format_read_test_all.o baseline_ncx_app.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) $(LIBS) -o $@ $^
	
h5_baseline_test_all: h5_baseline_test_all.o
	$(H5CC) $(CFLAGS) -o $@ $^ 

h5_%.o: h5_%.c
	$(H5CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

copy:
	cp ./lib_level_baseline_old.nc ./lib_level_baseline.nc
clean:
	rm -f $(OBJS)baseline_ex1 app_baseline_test_all pnc_consist_check lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all h5_baseline_test_all new_format_create_simple
base: baseline_test baseline_ex1 app_baseline_test_all pnc_consist_check
lib: $(LIB_PROGRAMS)
h5: h5_baseline_test_all
new_format: new_format_create_simple new_format_create_diff new_format_open new_format_create_large
check_format: new_format_create_simple new_format_open
	mpiexec -n 4 ./new_format_create_simple
	mpiexec -n 4 ./new_format_open
run_format:
	mpiexec -n 4 ./new_format_test_all
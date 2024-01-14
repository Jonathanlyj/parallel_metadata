CC = mpicc
CFLAGS = -O0 -g 
INCLUDES = -I/Users/liyoujia/Documents/CE499/pnetcdf/PnetCDF_install/include
LFLAGS = -L/Users/liyoujia/Documents/CE499/pnetcdf/PnetCDF_install/lib
LIBS = -lpnetcdf

SRCS = baseline_ex1.c baseline_ncx.c
OBJS = $(SRCS:.c=.o)

TARGET = baseline_ex1

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean :
	rm -f $(OBJS) $(TARGET)



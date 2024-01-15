CC = mpicc
CFLAGS = -O0 -g 
INCLUDES = -I$(PNETCDF_DIR)/include
LFLAGS = -L$(PNETCDF_DIR)/lib
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



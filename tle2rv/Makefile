FC = gfortran
CFLAGS = -c -O -fbounds-check
TARGET = get_iss_pv
OBJS = const.o ext.o model.o propagation.o io.o get_iss_pv.o
MODS = const.mod ext.mod model.mod propagation.mod io.mod

.SUFFIXES: .o .f95
.PHONY: all clean

all:	$(TARGET)

.f95.o:
	$(FC) $(CFLAGS) $<

$(TARGET):	$(OBJS)
	$(FC) -o $@ $(OBJS)

clean:
	@rm -f $(TARGET) $(OBJS) $(MODS)


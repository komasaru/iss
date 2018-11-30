FC = gfortran
CFLAGS = -c -O -fbounds-check
TARGET = iss_sgp4
OBJS = const.o ext.o model.o propagation.o io.o eph.o iss_sgp4.o
MODS = const.mod ext.mod model.mod propagation.mod io.mod eph.mod

.SUFFIXES: .o .f95
.PHONY: all clean

all:	$(TARGET)

.f95.o:
	$(FC) $(CFLAGS) $<

$(TARGET):	$(OBJS)
	$(FC) -o $@ $(OBJS)

clean:
	@rm -f $(TARGET) $(OBJS) $(MODS)


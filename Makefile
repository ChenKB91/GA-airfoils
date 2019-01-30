#Makefile for IBFS ADJOINT CODE MPI version
#(by Jeesoon Choi jeesoonchoi@gmail.com and Hsieh-Chen Tsai else731@gmail.com)

FFILES = parameters dts grid myfft
FFILES+= mmisc_bc mmisc_opr mmisc_motion mmisc
FFILES+= controls user write_responses
FFILES+= ibpm_pressure variables ibpm main

INFILES = ib
SHFILES = pbs_script_mpi
IBFILES = ${FFILES:=.o}

#COMPILE = /usr/local/mpich/mpich-3.2-gcc-7.1.0/bin/mpif90 #change
COMPILE = /usr/local/bin/gfortran
FLAGS1 =  -m64 -O3 -fcheck=bounds -I/usr/local/include
#FLAGS1 =  -m64 -O3 -heap-arrays -I/usr/local/fftw/fftw-3.3.6-pl2-gcc-7.1.0/include
#FLAGS1 =  -m64 -O3 -fcheck=bounds -I/usr/local/fftw/fftw-3.3.6-pl2-gcc-7.1.0/include #change
#FLAGS2 =  -L/usr/local/fftw/fftw-3.3.6-pl2-gcc-7.1.0/lib/ -lfftw3 -o  #change
FLAGS2 = -L/usr/local/lib/ -lfftw3 -o

BASE = /raid0/tflinois/BASE_adj

SRC = src
OBJ = obj
BIN = bin
OUT = output
IN  = input

DIR = ${SRC} ${OBJ} ${BIN} ${OUT} ${IN}
DIRS = ${DIR:=.dir}

IBFILES_SRC := $(foreach p, $(IBFILES), $(SRC)/$(p))
IBFILES_OBJ := $(foreach p, $(IBFILES), $(OBJ)/$(p))

#ib :  ${IBFILES_OBJ}
#	$(COMPILE) $(FLAGS1) $(IBFILES_OBJ) $(FLAGS2) $(BIN)/ib -J$(OBJ)
ib :  ${IBFILES_OBJ}
	$(COMPILE) $(FLAGS1) $(IBFILES_OBJ) $(FLAGS2) $(BIN)/ib

setup: ${DIRS}
	@for i in ${DIR} ; do \
               mkdir -p $${i} ;\
	done
	@for i in ${FFILES} ; do \
              echo "cp  -p ${BASE}/$${i}.f90 ${SRC}/$${i}.f90 " ;\
              cp -p ${BASE}/${SRC}/$${i}.f90 ${SRC}/ ;\
	done
	@for i in ${INFILES} ; do \
               echo "cp  -p ${BASE}/${IN}/$${i}.inp ${IN}/$${i}.inp " ;\
               cp  -p ${BASE}/${IN}/$${i}.inp ${IN}/$${i}.inp ;\
	done
	@for i in $(SHFILES) ; do \
               echo "cp  -p ${BASE}/$${i}.sh $${i}.sh" ;\
               cp  -p ${BASE}/$${i}.sh $${i}.sh ;\
        done

#realclean: ${DIRS} #
#	@for i in ${DIR} ; do \
#               rm -rf $${i} ;\
#	done

clean:
	rm -f $(BIN)/ib ${OBJ}/*.o *.mod

pyreset:
	python reset.py

cleandata:
	rm -f output/ib.chd output/responds/*.dat output/snapshots/*.var

#$(OBJ)/%.o: $(SRC)/%.f90
#	$(COMPILE) $(FLAGS1) -J$(@D) -o $@ -c $<
$(OBJ)/%.o: $(SRC)/%.f90
	$(COMPILE) $(FLAGS1) -o $@ -c $<

%.dir :
	$<

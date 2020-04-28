INCLUDE =
LDLIBSOPTIONS =
COPTIONS =

COPTIONS_DEBUG = -g -O0 -D_DEBUG -std=c++0x
COPTIONS_RELEASE = -O3 -std=c++0x

CC = mpicxx
LINK = mpicxx

SRC_DIR := src
BIN_DIR := bin

LINK_TARGET := libls

OBJS = triagGridTest.o

REBUILDABLES = $(SRC_DIR)/$(OBJS) $(BIN_DIR)/$(LINK_TARGET)

INCLUDE += \
	-I/home/petr/Libs/lsmlib-1.0.1/build_debug/include \
	-I/home/petr/Libs/eigen-3.1.2 \
	-I${PETSC_DIR}/include \
	-I${PETSC_DIR}/${PETSC_ARCH}/include

LDLIBSOPTIONS += \
	-L/home/petr/Libs/lsmlib-1.0.1/build_debug/lib \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib 

LDLIBSOPTIONS += -lboost_mpi -lboost_serialization -lpetsc -llapack -lblas -llsm_serial -llsm_toolbox -lm -lz

COPTIONS += ${INCLUDE}
COPTIONS += ${COPTIONS_RELEASE}

all: ${LINK_TARGET} 
	echo All done

# $@ expands to the rule's target, in this case "libls"
# $^ expands to the rule's dependencies, in this case petsctest.o
${LINK_TARGET}: ${SRC_DIR}/${OBJS}
	${CC} -o ${BIN_DIR}/$@ $^ ${LDLIBSOPTIONS}

${SRC_DIR}/%.o: ${SRC_DIR}/%.cpp
	${CC} -o $@ -c $< ${COPTIONS}

${SRC_DIR}/%.dep: ${SRC_DIR}/%.cpp
	${CC} ${COPTIONS} -M $< > $@ 
	include ${OBJS:.o=.dep}

clean:
	rm -f ${REBUILDABLES}

#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/AESCipher.o \
	${OBJECTDIR}/Approximation.o \
	${OBJECTDIR}/CombinatiorialGenerator.o \
	${OBJECTDIR}/ICipher.o \
	${OBJECTDIR}/ProgressMonitor.o \
	${OBJECTDIR}/base.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-lm -lcrypto -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp
CXXFLAGS=-lm -lcrypto -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation ${OBJECTFILES} ${LDLIBSOPTIONS} -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp

${OBJECTDIR}/AESCipher.o: AESCipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AESCipher.o AESCipher.cpp

${OBJECTDIR}/Approximation.o: Approximation.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Approximation.o Approximation.cpp

${OBJECTDIR}/CombinatiorialGenerator.o: CombinatiorialGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CombinatiorialGenerator.o CombinatiorialGenerator.cpp

${OBJECTDIR}/ICipher.o: ICipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ICipher.o ICipher.cpp

${OBJECTDIR}/ProgressMonitor.o: ProgressMonitor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ProgressMonitor.o ProgressMonitor.cpp

${OBJECTDIR}/base.o: base.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/base.o base.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Ifaugere -Ifaugere/int -Ifaugere/protocol -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

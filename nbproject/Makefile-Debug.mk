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
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/AESCipher.o \
	${OBJECTDIR}/AESCipherOpenSSL.o \
	${OBJECTDIR}/Approximation.o \
	${OBJECTDIR}/CombinatiorialGenerator.o \
	${OBJECTDIR}/CombinatorialIndexer.o \
	${OBJECTDIR}/FGbHelper.o \
	${OBJECTDIR}/ICipher.o \
	${OBJECTDIR}/Logger.o \
	${OBJECTDIR}/NTLUtils.o \
	${OBJECTDIR}/ProgressMonitor.o \
	${OBJECTDIR}/base.o \
	${OBJECTDIR}/ciphers/aes.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/sha3/Sha3Interface.o \
	${OBJECTDIR}/sha3/hash_functions/BMW/BMW_sha3.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakDuplex.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakF-1600-opt32.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakSponge.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak/Keccak_sha3.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -lboost_program_options -lcrypto -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp -lntl
CXXFLAGS=-m64 -lboost_program_options -lcrypto -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp -lntl

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-Lfaugere/x64 -lm -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf `pkg-config --libs libcrypto`  

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation ${OBJECTFILES} ${LDLIBSOPTIONS} -lgmp -lntl -lboost_program_options

${OBJECTDIR}/AESCipher.o: AESCipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AESCipher.o AESCipher.cpp

${OBJECTDIR}/AESCipherOpenSSL.o: AESCipherOpenSSL.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AESCipherOpenSSL.o AESCipherOpenSSL.cpp

${OBJECTDIR}/Approximation.o: Approximation.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Approximation.o Approximation.cpp

${OBJECTDIR}/CombinatiorialGenerator.o: CombinatiorialGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CombinatiorialGenerator.o CombinatiorialGenerator.cpp

${OBJECTDIR}/CombinatorialIndexer.o: CombinatorialIndexer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CombinatorialIndexer.o CombinatorialIndexer.cpp

${OBJECTDIR}/FGbHelper.o: FGbHelper.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FGbHelper.o FGbHelper.cpp

${OBJECTDIR}/ICipher.o: ICipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ICipher.o ICipher.cpp

${OBJECTDIR}/Logger.o: Logger.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Logger.o Logger.cpp

${OBJECTDIR}/NTLUtils.o: NTLUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/NTLUtils.o NTLUtils.cpp

${OBJECTDIR}/ProgressMonitor.o: ProgressMonitor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ProgressMonitor.o ProgressMonitor.cpp

${OBJECTDIR}/base.o: base.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/base.o base.cpp

${OBJECTDIR}/ciphers/aes.o: ciphers/aes.cpp 
	${MKDIR} -p ${OBJECTDIR}/ciphers
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ciphers/aes.o ciphers/aes.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/sha3/Sha3Interface.o: sha3/Sha3Interface.cpp 
	${MKDIR} -p ${OBJECTDIR}/sha3
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/Sha3Interface.o sha3/Sha3Interface.cpp

${OBJECTDIR}/sha3/hash_functions/BMW/BMW_sha3.o: sha3/hash_functions/BMW/BMW_sha3.cpp 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/BMW
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/BMW/BMW_sha3.o sha3/hash_functions/BMW/BMW_sha3.cpp

${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakDuplex.o: sha3/hash_functions/Keccak/KeccakDuplex.c 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak
	${RM} "$@.d"
	$(COMPILE.c) -g `pkg-config --cflags libcrypto`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakDuplex.o sha3/hash_functions/Keccak/KeccakDuplex.c

${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakF-1600-opt32.o: sha3/hash_functions/Keccak/KeccakF-1600-opt32.c 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak
	${RM} "$@.d"
	$(COMPILE.c) -g `pkg-config --cflags libcrypto`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakF-1600-opt32.o sha3/hash_functions/Keccak/KeccakF-1600-opt32.c

${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakSponge.o: sha3/hash_functions/Keccak/KeccakSponge.c 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak
	${RM} "$@.d"
	$(COMPILE.c) -g `pkg-config --cflags libcrypto`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak/KeccakSponge.o sha3/hash_functions/Keccak/KeccakSponge.c

${OBJECTDIR}/sha3/hash_functions/Keccak/Keccak_sha3.o: sha3/hash_functions/Keccak/Keccak_sha3.cpp 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak/Keccak_sha3.o sha3/hash_functions/Keccak/Keccak_sha3.cpp

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

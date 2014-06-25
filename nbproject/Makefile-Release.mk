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
	${OBJECTDIR}/AESCipherOpenSSL.o \
	${OBJECTDIR}/Approximation.o \
	${OBJECTDIR}/CombinatiorialGenerator.o \
	${OBJECTDIR}/CombinatorialIndexer.o \
	${OBJECTDIR}/FGbHelper.o \
	${OBJECTDIR}/ICipher.o \
	${OBJECTDIR}/Keccak2.o \
	${OBJECTDIR}/KeccakFull.o \
	${OBJECTDIR}/KeccakOptAsm7r.o \
	${OBJECTDIR}/Logger.o \
	${OBJECTDIR}/MultiCombinatorialGenerator.o \
	${OBJECTDIR}/NTLUtils.o \
	${OBJECTDIR}/ProgressMonitor.o \
	${OBJECTDIR}/base.o \
	${OBJECTDIR}/ciphers/aes.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/sha3/Sha3Interface.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakF-1600-x86-64-asm.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakSponge.o \
	${OBJECTDIR}/sha3/hash_functions/Keccak64_r7/KeccakF-1600-x86-64-gas.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -lboost_program_options -lcrypto -dynamic -ldl -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp -lntl -fomit-frame-pointer -O3 -g0 -march=native -mtune=native -m64
CXXFLAGS=-m64 -lboost_program_options -lcrypto -dynamic -ldl -L faugere/x64 -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp -lm -fopenmp -lntl -fomit-frame-pointer -O3 -g0 -march=native -mtune=native -m64

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
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/highorderapproximation ${OBJECTFILES} ${LDLIBSOPTIONS} -lgmp -lntl -lboost_program_options -lboost_iostreams -lboost_serialization -ldl -lz

${OBJECTDIR}/AESCipher.o: AESCipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AESCipher.o AESCipher.cpp

${OBJECTDIR}/AESCipherOpenSSL.o: AESCipherOpenSSL.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AESCipherOpenSSL.o AESCipherOpenSSL.cpp

${OBJECTDIR}/Approximation.o: Approximation.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Approximation.o Approximation.cpp

${OBJECTDIR}/CombinatiorialGenerator.o: CombinatiorialGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CombinatiorialGenerator.o CombinatiorialGenerator.cpp

${OBJECTDIR}/CombinatorialIndexer.o: CombinatorialIndexer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CombinatorialIndexer.o CombinatorialIndexer.cpp

${OBJECTDIR}/FGbHelper.o: FGbHelper.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FGbHelper.o FGbHelper.cpp

${OBJECTDIR}/ICipher.o: ICipher.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ICipher.o ICipher.cpp

${OBJECTDIR}/Keccak2.o: Keccak2.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Keccak2.o Keccak2.cpp

${OBJECTDIR}/KeccakFull.o: KeccakFull.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/KeccakFull.o KeccakFull.cpp

${OBJECTDIR}/KeccakOptAsm7r.o: KeccakOptAsm7r.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/KeccakOptAsm7r.o KeccakOptAsm7r.cpp

${OBJECTDIR}/Logger.o: Logger.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Logger.o Logger.cpp

${OBJECTDIR}/MultiCombinatorialGenerator.o: MultiCombinatorialGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MultiCombinatorialGenerator.o MultiCombinatorialGenerator.cpp

${OBJECTDIR}/NTLUtils.o: NTLUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/NTLUtils.o NTLUtils.cpp

${OBJECTDIR}/ProgressMonitor.o: ProgressMonitor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ProgressMonitor.o ProgressMonitor.cpp

${OBJECTDIR}/base.o: base.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/base.o base.cpp

${OBJECTDIR}/ciphers/aes.o: ciphers/aes.cpp 
	${MKDIR} -p ${OBJECTDIR}/ciphers
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ciphers/aes.o ciphers/aes.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/sha3/Sha3Interface.o: sha3/Sha3Interface.cpp 
	${MKDIR} -p ${OBJECTDIR}/sha3
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -I. -Ifaugere -Ifaugere/int -Ifaugere/protocol `pkg-config --cflags libcrypto` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/Sha3Interface.o sha3/Sha3Interface.cpp

${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakF-1600-x86-64-asm.o: sha3/hash_functions/Keccak64_common/KeccakF-1600-x86-64-asm.c 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak64_common
	${RM} "$@.d"
	$(COMPILE.c) -O2 `pkg-config --cflags libcrypto`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakF-1600-x86-64-asm.o sha3/hash_functions/Keccak64_common/KeccakF-1600-x86-64-asm.c

${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakSponge.o: sha3/hash_functions/Keccak64_common/KeccakSponge.c 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak64_common
	${RM} "$@.d"
	$(COMPILE.c) -O2 `pkg-config --cflags libcrypto`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sha3/hash_functions/Keccak64_common/KeccakSponge.o sha3/hash_functions/Keccak64_common/KeccakSponge.c

${OBJECTDIR}/sha3/hash_functions/Keccak64_r7/KeccakF-1600-x86-64-gas.o: sha3/hash_functions/Keccak64_r7/KeccakF-1600-x86-64-gas.s 
	${MKDIR} -p ${OBJECTDIR}/sha3/hash_functions/Keccak64_r7
	$(AS) $(ASFLAGS) -o ${OBJECTDIR}/sha3/hash_functions/Keccak64_r7/KeccakF-1600-x86-64-gas.o sha3/hash_functions/Keccak64_r7/KeccakF-1600-x86-64-gas.s

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

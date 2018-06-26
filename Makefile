
ifdef COMSPEC
  CC = icl
  SRC = src\*.c
  LIB = -Qopenmp -Qmkl
  OPTIONS = -Qstd=c99 -MT
  MIC_ENABLED = -Qmic

  INTEL64LIB = "E:\Intel\Composer XE 2013\compiler\lib\intel64"
  VCLIB = "E:\VisualStudio\VC\lib\amd64"
  WINSDK = "C:\Program Files (x86)\Windows Kits\8.0\Lib\win8\um\x64"

  LIB_PATH = -link -libpath:$(INTEL64LIB)  -link -libpath:$(VCLIB)  -link -libpath:$(WINSDK)
else
  CC = icc
  SRC = $(wildcard src/*.c)
  LIB = -lm -qopenmp -mkl -D_GNU_SOURCE
  OPTIONS = -Wall -O3
  MIC_ENABLED = -mmic
  LIB_PATH =
endif

BUILD_TIME := $(shell date)
GIT_TAG := $(shell git describe --abbrev=0 --tags)
GIT_BRANCH := $(shell git branch | grep \* | cut -d ' ' -f2-)
GIT_COMMIT := $(shell git rev-parse HEAD)
GIT_COMMIT_TIME := $(shell git show -s --format=%ci)
GIT_DESCR := $(shell git describe --dirty --always --tags)
FULL_VERSION = -DVERSION="\"Build date: $(BUILD_TIME)\nVersion: $(GIT_TAG)\nGit branch: $(GIT_BRANCH)\nCommit: $(GIT_COMMIT)\nCommit time: $(GIT_COMMIT_TIME)\nGit description: $(GIT_DESCR)\""

default : MCsquare_linux MCsquare_linux_NoArch

MCsquare_linux : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) $(FULL_VERSION) -static -m64 -march=corei7-avx -o MCsquare_linux

MCsquare_linux_NoArch : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) $(FULL_VERSION) -static -m64 -o MCsquare_linux_NoArch

MCsquare_linux_RC :	$(SRC)
	make MCsquare_linux
	mv ./MCsquare_linux ./MCsquare_linux_RC

MCsquare_mac :	$(SRC)
	$(CC) $(SRC) -mkl -lm -qopenmp $(OPTIONS)  -m64 -march=corei7-avx -D_FORTIFY_SOURCE=1 -o MCsquare_mac

MCsquare_mac_NoArch :	$(SRC)
	$(CC) $(SRC) -mkl -lm -qopenmp $(OPTIONS)  -m64 -D_FORTIFY_SOURCE=1 -o MCsquare_mac_NoArch

MC2_gcc_UMCG : $(SRC)
	gcc $(SRC) -fcilkplus -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lm -ldl $(OPTIONS) $(LIB_PATH) -static -m64 -march=corei7-avx  -I/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/include/ -L/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64_lin -o MC2_gcc





clean:
	rm -f *~
	rm -f src/*~
	rm -f src/include/*~
	rm -f src/include/lib/*~

lines:
	make clean
	find src -type f -exec cat {} + | wc -l

shared : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -m64 -march=corei7-avx -o shared

debug : $(SRC)
	$(CC) $(SRC) $(LIB) $(LIB_PATH) -static -g -debug inline-debug-info -o debug

profile : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -p -g -o profile

MC2_gcc : $(SRC)
	gcc $(SRC) -fcilkplus -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_vml_def -lm $(OPTIONS) $(LIB_PATH) -m64 -march=corei7-avx -o MC2_gcc

debug_gcc : $(SRC)
	gcc $(SRC) -fcilkplus -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_vml_def -lm $(OPTIONS) $(LIB_PATH) -g -ggdb -rdynamic -o debug_gcc



debug_fp : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -static -fp-trap-all=common -traceback -g -o debug_fp

mic : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(MIC_ENABLED) $(LIB_PATH) -static -o mic

mic_shared : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(MIC_ENABLED) $(LIB_PATH) -static -o mic_shared

report : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -static -m64 -march=corei7-avx -vec-report=3 -O3 -o report

mic_report : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(MIC_ENABLED) $(LIB_PATH) -static -vec-report=3 -o mic_report

novec : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(MIC_ENABLED) $(LIB_PATH) -static -no-vec -o novec

novec_cpu : $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -static -no-vec -o novec_cpu

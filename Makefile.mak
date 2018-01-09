
CC = icl
SRC = src\*.c
LIB = -Qopenmp -Qmkl
OPTIONS = -Qstd=c99 -MT
MIC_ENABLED = -Qmic

INTEL64LIB = "E:\Intel\Composer XE 2013\compiler\lib\intel64"
VCLIB = "E:\VisualStudio\VC\lib\amd64"
WINSDK = "C:\Program Files (x86)\Windows Kits\8.0\Lib\win8\um\x64"

LIB_PATH = -link -libpath:$(INTEL64LIB)  -link -libpath:$(VCLIB)  -link -libpath:$(WINSDK)


MCsquare_win.exe: $(SRC)
	$(CC) $(SRC) $(LIB) $(OPTIONS) $(LIB_PATH) -link -out:MCsquare_win.exe

@echo off

del MCsquare_win*

echo CC = icl > Makfile.mak
echo SRC = src\*.c >> Makfile.mak
echo LIB = -Qopenmp -Qmkl >> Makfile.mak
echo OPTIONS = -Qstd=c99 -MT -Qdiag-disable:10441 /O3 >> Makfile.mak
echo MIC_ENABLED = -Qmic >> Makfile.mak

echo INTEL64LIB = "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win" >> Makfile.mak
echo MKLIB = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\intel64" >> Makfile.mak
echo VCLIB = "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.32.31326\lib\x64" >> Makfile.mak
echo WINSDK1 = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.19041.0\ucrt\x64" >> Makfile.mak
echo WINSDK2 = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.19041.0\um\x64" >> Makfile.mak
echo LIB_PATH = -link -libpath:$(INTEL64LIB)  -link -libpath:$(MKLIB)  -link -libpath:$(VCLIB)  -link -libpath:$(WINSDK1)  -link -libpath:$(WINSDK2) >> Makfile.mak


for /f "delims=" %%i in ('date /t') do set BUILD_TIME=%%i
for /f "delims=" %%i in ('git describe --abbrev^=0 --tags') do set GIT_TAG=%%i
for /f "delims=" %%i in ('git branch ^| findstr *') do set GIT_BRANCH=%%i
for /f "delims=" %%i in ('git rev-parse HEAD') do set GIT_COMMIT=%%i
for /f "delims=" %%i in ('git show -s --format^=%%ci') do set GIT_COMMIT_TIME=%%i
for /f "delims=" %%i in ('git describe --dirty --always --tags') do set GIT_DESCR=%%i

echo FULL_VERSION = /DVERSION="\"Build date: %BUILD_TIME%\nVersion: %GIT_TAG%\nGit branch: %GIT_BRANCH%\nCommit: %GIT_COMMIT%\nCommit time: %GIT_COMMIT_TIME%\nGit description: %GIT_DESCR%\"" >> Makfile.mak

set "TAB=	"

echo default : MCsquare_win.exe >> Makfile.mak
echo: >> Makfile.mak
echo all : MCsquare_win.exe MCsquare_win_sse4.exe MCsquare_win_avx.exe MCsquare_win_avx2.exe MCsquare_win_avx512.exe >> Makfile.mak
echo: >> Makfile.mak

echo MCsquare_win.exe: $(SRC) >> Makfile.mak
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win.exe >> Makfile.mak
echo: >> Makfile.mak

echo MCsquare_win_sse4.exe: $(SRC) >> Makfile.mak
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) /QxSSE4.2 $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win_sse4.exe >> Makfile.mak
echo: >> Makfile.mak

echo MCsquare_win_avx.exe: $(SRC) >> Makfile.mak
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) /QxAVX $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win_avx.exe >> Makfile.mak
echo: >> Makfile.mak

echo MCsquare_win_avx2.exe: $(SRC) >> Makfile.mak
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) /QxCORE-AVX2 $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win_avx2.exe >> Makfile.mak
echo: >> Makfile.mak

echo MCsquare_win_avx512.exe: $(SRC) >> Makfile.mak
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) /QxCORE-AVX512 $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win_avx512.exe >> Makfile.mak
echo: >> Makfile.mak


nmake -f Makfile.mak all

del *.obj

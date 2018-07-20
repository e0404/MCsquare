@echo off


echo CC = icl > Makfile.mak
echo SRC = src\*.c >> Makfile.mak
echo LIB = -Qopenmp -Qmkl >> Makfile.mak
echo OPTIONS = -Qstd=c99 -MT >> Makfile.mak
echo MIC_ENABLED = -Qmic >> Makfile.mak

echo INTEL64LIB = "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\compiler\lib\intel64" >> Makfile.mak
echo VCLIB = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\lib\x64" >> Makfile.mak
echo WINSDK1 = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\ucrt\x64" >> Makfile.mak
echo WINSDK2 = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\um\x64" >> Makfile.mak
echo LIB_PATH = -link -libpath:$(INTEL64LIB)  -link -libpath:$(VCLIB)  -link -libpath:$(WINSDK1)  -link -libpath:$(WINSDK2) >> Makfile.mak


for /f "delims=" %%i in ('date /t') do set BUILD_TIME=%%i
for /f "delims=" %%i in ('git describe --abbrev^=0 --tags') do set GIT_TAG=%%i
for /f "delims=" %%i in ('git branch ^| findstr *') do set GIT_BRANCH=%%i
for /f "delims=" %%i in ('git rev-parse HEAD') do set GIT_COMMIT=%%i
for /f "delims=" %%i in ('git show -s --format^=%%ci') do set GIT_COMMIT_TIME=%%i
for /f "delims=" %%i in ('git describe --dirty --always --tags') do set GIT_DESCR=%%i

echo FULL_VERSION = /DVERSION="\"Build date: %BUILD_TIME%\nVersion: %GIT_TAG%\nGit branch: %GIT_BRANCH%\nCommit: %GIT_COMMIT%\nCommit time: %GIT_COMMIT_TIME%\nGit description: %GIT_DESCR%\"" >> Makfile.mak

echo MCsquare_win.exe: $(SRC) >> Makfile.mak

set "TAB=	"
echo %TAB% $(CC) $(SRC) $(LIB) $(OPTIONS) $(FULL_VERSION) $(LIB_PATH) -link -out:MCsquare_win.exe >> Makfile.mak


nmake -f Makfile.mak
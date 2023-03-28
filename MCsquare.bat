@echo off

SET This_Dir=%~dp0
C:\Windows\System32\setx MCsquare_Materials_Dir %This_Dir%Materials >NUL

echo MCsquare is running ...

call %This_Dir%MCsquare_win_avx512.exe | findstr /i "processor support Intel(R)" >NUL || exit 
call %This_Dir%MCsquare_win_avx2.exe | findstr /i "processor support Intel(R)" >NUL || exit
call %This_Dir%MCsquare_win_avx.exe | findstr /i "processor support Intel(R)" >NUL || exit
call %This_Dir%MCsquare_win_sse4.exe | findstr /i "processor support Intel(R)" >NUL || exit
call %This_Dir%MCsquare_win.exe 


ECHO OFF

REM extract the input.zip folder
.\7za.exe x -o.\ input%1.zip

REM move into container FOLDER
@for /f "tokens=* delims= " %%a in ('@dir /ad /b') do @set FOLDER=%%a

REM extract IM.zip into container FOLDER
.\7za.exe x -o.\%FOLDER%\ Research.zip

cd .\%FOLDER%

.\shell.exe

cd ..

REM create the output.zip file. customize the file list for each program
.\7za.exe a -tzip %FOLDER%.zip .\%FOLDER%\output


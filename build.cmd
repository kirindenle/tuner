@echo off

IF NOT EXIST build mkdir build
IF NOT EXIST build\SDL2.dll copy external\SDL2-2.0.12\lib\x64\SDL2.dll build\SDL2.dll
IF NOT EXIST build\SDL2_ttf.dll copy external\SDL2_ttf-2.0.15\lib\x64\SDL2_ttf.dll build\SDL2_ttf.dll
IF NOT EXIST build\libfreetype-6.dll copy external\SDL2_ttf-2.0.15\lib\x64\libfreetype-6.dll build\libfreetype-6.dll

SETLOCAL
if "%~1"=="" (
    set Name="game1"
) else (
    set Name="%~1"
)
set includeSDL=external\SDL2-2.0.12\include
set includeSDL_ttf=external\SDL2_ttf-2.0.15\include
set includeCSTDLIB="C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\ucrt"
set includeCRuntime="C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Tools\MSVC\14.16.27023\include"
set includeMyLib="C:\Users\22599\projects\_mylib"
set SDLlibs=SDL2main.lib SDL2.lib Shell32.lib SDL2_ttf.lib
set SDLlibpath=external\SDL2-2.0.12\lib\x64
set SDL_ttf_libpath=external\SDL2_ttf-2.0.15\lib\x64
set SourceCode=src/%Name%.cpp
set Executable=build/%Name%.exe
set Obj=build/%Name%.obj
set Pdb=build/%Name%.pdb

call cl /nologo /EHsc %PreprocessorDefines% /I %includeSDL% /I %includeSDL_ttf% /I %includeCSTDLIB% /I %includeCRuntime% /I %includeMyLib% ^
    /Z7 ^
    /W4 /wd4189 ^
    /std:c++17 ^
    /Fe:%Executable% /Fo:%Obj% ^
    /Tp %SourceCode% ^
    /link /LIBPATH:%SDLlibpath% /LIBPATH:%SDL_ttf_libpath% %SDLlibs% ^
    /PDB:%Pdb% ^
    /SUBSYSTEM:WINDOWS
ENDLOCAL

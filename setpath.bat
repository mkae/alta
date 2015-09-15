:: BAT command to set ALTA environment for ease of use.
:: Invoke this script on Windows to have access to ALTA's
:: commands in the cmd.exe
::
::    setpath.bat
::
:: ALTA command line programs and plugins will be available
:: to the shell.
::
:: This script can be launch using Administrator rights. In
:: this case, the different PATHS will be permanent.
::
net session >nul 2>&1
if %ERRORLEVEL% equ 0 (
	setx ALTA_DIR "%~dp0sources" /M
	setx ALTA_LIB "%~dp0sources\build" /M
	setx PATH "%PATH%;%~dp0sources\build" /M
	setx PYTHONPATH "%PYTHONPATH%;%~dp0sources\build" /M
) else (
	set ALTA_DIR=%~dp0sources\
	set ALTA_LIB=%~dp0sources\build
	set PATH=%PATH%;%~dp0sources\build
	set PYTHONPATH=%PYTHONPATH%;%~dp0sources\build
)
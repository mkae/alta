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


:: It can happens that the PATH variable contain directory
:: with parenthesis in them. In that case, we cannot use
:: a parenthesis for the IF ELSE statement. I use goto instead.
::
:: Source:
::   http://www.blinnov.com/en/2010/06/04/microsoft-was-unexpected-at-this-time/
net session >nul 2>&1
IF %ERRORLEVEL% EQU 0 GOTO superuser

set ALTA_DIR=%~dp0sources
set ALTA_PLUGIN_PATH=%~dp0build\plugins
set PATH=%PATH%;%~dp0build\softs
set PYTHONPATH=%PYTHONPATH%;%~dp0build\python
REM set ALTA_DIR=%CD%\sources
REM set ALTA_PLUGIN_PATH=%CD%\build\plugins
REM set PATH=%PATH%;%CD%\build\softs
REM set PYTHONPATH=%PYTHONPATH%;%CD%\build\python
GOTO:eof

:superuser 
setx ALTA_DIR "%~dp0sources" /M
setx ALTA_PLUGIN_PATH "%~dp0build\plugins" /M
setx PATH "%PATH%;%~dp0build\softs" /M
setx PYTHONPATH "%PYTHONPATH%;%~dp0build\python" /M
GOTO:eof

:: BAT command to set ALTA environment for ease of use.
:: Invoke this script on Windows to have access to ALTA's
:: commands in the cmd.exe
::
::    setpath.bat
::
:: ALTA command line programs and plugins will be available
:: to the shell.
::
:: NOTE: This needs to be launched at the root of ALTA. If
:: not, the programs will not be accessible.
::
set ALTA_DIR=%cd%\sources\
set ALTA_LIB=%cd%\sources\build
set PATH=%PATH%;%cd%\sources\build
set PYTHONPATH=%PYTHONPATH%;%cd%\sources\build

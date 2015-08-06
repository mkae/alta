!include EnvVarUpdate.nsh

# Name of the installer
Outfile "ALTA.exe"

InstallDir "$DESKTOP\ALTA"

!define ALTADIR "${__FILEDIR__}\..\.."

# Request application privileges for Windows Vista
RequestExecutionLevel user

# Section start
Section

SetOutPath $INSTDIR\bin
File "${ALTADIR}\sources\build\*.exe"

SetOutPath $INSTDIR\lib
File "${ALTADIR}\sources\build\core.lib"

SetOutPath $INSTDIR\plugins
File "${ALTADIR}\sources\build\nonlinear*.dll"
File "${ALTADIR}\sources\build\rational*.dll"
File "${ALTADIR}\sources\build\data*.dll"

SetOutPath $INSTDIR\python
File "${ALTADIR}\sources\build\alta.dll"

# Update the ENVIROMNENT
WriteRegStr HKCU "Environment" "ALTA_DIR"   '$INSTDIR'
WriteRegStr HKCU "Environment" "ALTA_LIB"   '$INSTDIR\plugins'
${EnvVarUpdate} $0 "PATH"     "A" "HKCU" '$INSTDIR\bin'
${EnvVarUpdate} $1 "PYTHONPATH" "A" "HKCU" '$INSTDIR\python'

#End Section
SectionEnd

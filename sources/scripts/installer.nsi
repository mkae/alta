!include EnvVarUpdate.nsh
!include "MUI2.nsh"
!include "winmessages.nsh"


# Name of the installer
Outfile "ALTA.exe"

InstallDir "$DESKTOP\ALTA"

!define ALTADIR "${__FILEDIR__}\..\.."

# Request application privileges for Windows Vista
RequestExecutionLevel user


# Pages
  !insertmacro MUI_PAGE_WELCOME
  !define MUI_WELCOMEPAGE_TITLE "Welcome to ALTA Installer"

  !insertmacro MUI_PAGE_LICENSE "..\..\LICENSE.txt"
#	!insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES




  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES

;Interface Settings

  !define MUI_ABORTWARNING

  !insertmacro MUI_LANGUAGE "English"


# Section start
Section "ALTA" SecMain

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



  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

#End Section
SectionEnd

;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;FILES

  Delete "$INSTDIR\Uninstall.exe"

  ;Remove the directory and all its content!
  RMDir /r "$INSTDIR"


  ;Remove ALTA_LIB AND ALTA_DIR env. variable
  DeleteRegKey HKCU "ALTA_DIR"
  DeleteRegKey HKCU "ALTA_LIB"


  ${un.EnvVarUpdate} $0 "PATH" "R" "HKCU" "$INSTDIR"
  ${un.EnvVarUpdate} $0 "PYTHONPATH" "R" "HKCU" "$INSTDIR\python"


  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

SectionEnd

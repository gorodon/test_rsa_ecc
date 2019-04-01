set Configuration=Release
set VTOOLSET=14.0

:: ----------------------------------------------------------------------

set BPlatform=Win32

set Out_Dir=.\%Configuration%

set MSBUILD="C:\Program Files (x86)\MSBuild\%VTOOLSET%\Bin\MSBuild.exe"

:: Build start

%MSBUILD% ./test_rsa_raw.vcxproj /p:Configuration=%Configuration%;OutDir=%Out_Dir%\ /p:Platform="%BPlatform%" /toolsversion:14.0 /t:rebuild

:: ----------------------------------------------------------------------

set BPlatform=x64

set Out_Dir=.\%Configuration%_%BPlatform%

set MSBUILD="C:\Program Files (x86)\MSBuild\%VTOOLSET%\Bin\amd64\MSBuild.exe"

:: Build start

%MSBUILD% ./test_rsa_raw.vcxproj /p:Configuration=%Configuration%;OutDir=%Out_Dir%\ /p:Platform="%BPlatform%" /toolsversion:14.0 /t:rebuild
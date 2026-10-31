$ErrorActionPreference = "Stop"
# FLAC3D 6.0 GQ/H-MRDF UDM build. Flags follow the GQ/H v6 port (SDK headers are 2016-era:
# /std:c++14, no /permissive-, forced iostream/string/vector includes, /W3).
# The SDK example project is VS2010 format (ToolsVersion 4.0, 10.0.40219.1); the install also
# ships a VS2015 add-in. MSVC v142 (VS2019) is binary compatible with the VS2015 C++ ABI.
$root = Split-Path -Parent $PSScriptRoot
# FLAC3D 6 SDK: set FLAC3D6_PLUGINFILES, otherwise the default install location is used.
$sdk = if ($env:FLAC3D6_PLUGINFILES) { $env:FLAC3D6_PLUGINFILES } else { Join-Path $env:ProgramFiles "Itasca\Flac3d600\pluginfiles" }
if (-not (Test-Path -LiteralPath $sdk)) { throw "FLAC3D 6 pluginfiles not found: $sdk (set FLAC3D6_PLUGINFILES)" }
$pf86 = ${env:ProgramFiles(x86)}
$vcvars = Join-Path $pf86 "Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path -LiteralPath $vcvars)) {
    $vswhere = Join-Path $pf86 "Microsoft Visual Studio\Installer\vswhere.exe"
    $vsroot = & $vswhere -version "[16.0,17.0)" -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    $vcvars = Join-Path $vsroot "VC\Auxiliary\Build\vcvars64.bat"
}
if (-not (Test-Path -LiteralPath $vcvars)) { throw "vcvars64.bat not found" }
$src = Join-Path $root "src_v6\modelgqh_mrdf_v6.cpp"
$dll = Join-Path $PSScriptRoot "cmodelgqh-mrdf006.dll"
$lib = Join-Path $PSScriptRoot "cmodelgqh-mrdf006.lib"
$inc = @("/I`"$sdk`"", "/I`"$sdk\cmodels`"", "/I`"$sdk\cmodels\src`"", "/I`"$sdk\interface`"", "/I`"$sdk\interface\base\src`"") -join " "
$cmd = "call `"$vcvars`" && cd /d `"$PSScriptRoot`" && cl.exe /nologo /std:c++14 /EHsc /LD /O2 /fp:precise /W3 /MD /DWIN32 /DNDEBUG /DDIM=3 /FI`"iostream`" /FI`"string`" /FI`"vector`" $inc /Fe:`"$dll`" `"$src`" /link /LIBPATH:`"$sdk\lib\exe64`" conmodel006_64.lib base005_64.lib /OUT:`"$dll`" /IMPLIB:`"$lib`""
Write-Output $cmd
cmd.exe /d /s /c $cmd
if ($LASTEXITCODE -ne 0) { throw "build failed: $LASTEXITCODE" }
Write-Output $dll

$ErrorActionPreference = "Stop"

$target = Split-Path -Parent $PSScriptRoot
$sdk = if ($env:FLAC3D9_PLUGINFILES) { $env:FLAC3D9_PLUGINFILES } else { "C:\ProgramData\Itasca\ItascaSoftware900\PluginFiles" }
# MSVC'yi bul: once 2019 BuildTools'un bilinen yolu, yoksa vswhere ile en yeni
# kurulum (2022 Build Tools dahil). Toolset v142/v143 fark etmez; ikisi de
# VS2015+ C++ ABI'siyle uyumludur.
$vcvars = "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path -LiteralPath $vcvars)) {
    $vcvars = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
}
if (-not (Test-Path -LiteralPath $vcvars)) {
    $vswhere = "C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path -LiteralPath $vswhere) {
        $vsroot = & $vswhere -latest -products * `
            -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
            -property installationPath
        if ($vsroot) {
            $vcvars = Join-Path $vsroot "VC\Auxiliary\Build\vcvars64.bat"
        }
    }
}
$source = Join-Path $target "src\modelgqh_mrdf.cpp"
$refcurveSource = Join-Path $target "src\refcurve.cpp"
$dll = Join-Path $PSScriptRoot "cmodelgqh-mrdf009.dll"
$importLib = Join-Path $PSScriptRoot "cmodelgqh-mrdf009.lib"
$pdb = Join-Path $PSScriptRoot "cmodelgqh-mrdf009.pdb"

if (-not (Test-Path -LiteralPath $vcvars)) {
    throw "MSVC 2019 vcvars64.bat not found: $vcvars"
}
if (-not (Test-Path -LiteralPath (Join-Path $sdk "lib\exe64\conmodel009.lib"))) {
    throw "FLAC3D v9 conmodel009.lib not found under $sdk."
}

$includeFlags = @(
    "/I`"$sdk`"",
    "/I`"$(Join-Path $sdk 'interface')`"",
    "/I`"$(Join-Path $sdk 'models\src')`""
) -join " "
$libDir = Join-Path $sdk "lib\exe64"
$command = "call `"$vcvars`" && cl.exe /nologo /std:c++17 /EHsc /LD /O2 /fp:precise /W4 /WX /permissive- /DWIN32 /DNDEBUG /DDIM=3 $includeFlags /Fe:`"$dll`" /Fd:`"$pdb`" `"$source`" `"$refcurveSource`" /link /DEBUG /LIBPATH:`"$libDir`" conmodel009.lib base009.lib /OUT:`"$dll`" /IMPLIB:`"$importLib`""
Write-Output $command
cmd.exe /d /s /c $command
if ($LASTEXITCODE -ne 0) {
    throw "FLAC3D v9 UDM compilation failed with exit code $LASTEXITCODE."
}

Write-Output $dll

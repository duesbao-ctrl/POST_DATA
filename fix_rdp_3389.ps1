[CmdletBinding()]
param()

function Test-IsAdministrator {
    $currentIdentity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($currentIdentity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Ensure-Administrator {
    if (Test-IsAdministrator) {
        return
    }

    Write-Host 'Requesting administrator privileges...' -ForegroundColor Yellow
    $self = $PSCommandPath
    if ([string]::IsNullOrWhiteSpace($self)) {
        throw 'Unable to determine script path for elevation.'
    }

    $escaped = '"' + $self.Replace('"', '""') + '"'
    Start-Process -FilePath 'powershell.exe' -Verb RunAs -ArgumentList @(
        '-ExecutionPolicy', 'Bypass',
        '-File', $escaped
    )
    exit
}

Ensure-Administrator

$ErrorActionPreference = 'Stop'

Write-Host 'Enabling Remote Desktop...' -ForegroundColor Cyan
Set-ItemProperty -Path 'HKLM:\SYSTEM\CurrentControlSet\Control\Terminal Server' -Name fDenyTSConnections -Value 0

Write-Host 'Ensuring RDP port is 3389...' -ForegroundColor Cyan
Set-ItemProperty -Path 'HKLM:\SYSTEM\CurrentControlSet\Control\Terminal Server\WinStations\RDP-Tcp' -Name PortNumber -Value 3389

Write-Host 'Enabling Remote Desktop firewall rules...' -ForegroundColor Cyan
Enable-NetFirewallRule -DisplayGroup 'Remote Desktop' | Out-Null

Write-Host 'Configuring and starting Remote Desktop Services...' -ForegroundColor Cyan
Set-Service -Name TermService -StartupType Automatic
Start-Service -Name TermService

try {
    Set-Service -Name UmRdpService -StartupType Manual
    Start-Service -Name UmRdpService -ErrorAction SilentlyContinue
} catch {
    Write-Host 'UmRdpService could not be started explicitly; continuing.' -ForegroundColor DarkYellow
}

Start-Sleep -Seconds 2

$listener = Get-NetTCPConnection -State Listen -ErrorAction SilentlyContinue |
    Where-Object { $_.LocalPort -eq 3389 } |
    Select-Object -First 1

$termService = Get-Service -Name TermService
$rdpEnabled = (Get-ItemProperty -Path 'HKLM:\SYSTEM\CurrentControlSet\Control\Terminal Server').fDenyTSConnections -eq 0

Write-Host ''
Write-Host 'RDP status:' -ForegroundColor Green
Write-Host ("  Remote Desktop enabled : {0}" -f $rdpEnabled)
Write-Host ("  TermService status     : {0}" -f $termService.Status)

if ($null -ne $listener) {
    Write-Host ("  Listening endpoint     : {0}:{1}" -f $listener.LocalAddress, $listener.LocalPort)
    Write-Host ''
    Write-Host 'RDP listener is ready. natfrp can forward 127.0.0.1:3389 now.' -ForegroundColor Green
} else {
    Write-Host '  Listening endpoint     : NOT LISTENING' -ForegroundColor Red
    Write-Host ''
    Write-Host 'RDP is still not listening on 3389. A reboot or additional policy check may be required.' -ForegroundColor Yellow
    exit 1
}

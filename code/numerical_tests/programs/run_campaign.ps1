# run_campaign.ps1 -- launches the campaign processes of the first-draft numerical test.
# Each process: python tbx_run.py OUT FAMILY K SEED NPACK DELTA THETAS NY TLIMIT  (internal limit TLIMIT s,
# 380 MB working-set guard inside the flow tracer).  At most $MaxPar processes at a time; a process is started
# only if more than 800 MB of physical memory are free (otherwise wait 20 s, up to 15 times, then skip);
# a process still running after 570 s is stopped by its PID.  Everything is logged to out/resource_log.txt.
# Campaign A (2026-10-07) used the defaults below with -Tag A -NPack 40 and thetas chain,a+ (random: chain,0.02,a+).
param([string]$Tag = "A", [int]$MaxPar = 3, [int]$NPack = 40, [int]$NY = 24, [int]$TLimit = 420,
      [string]$Fams = "random,rows,column,phase,lshape,rotgrid,stair,converge", [string]$Ks = "14,16",
      [string]$Deltas = "1/20,1/10", [string]$Thetas = "", [int]$SeedBase = 1000)
$py = "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe"
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force out | Out-Null
$log = "out/resource_log.txt"
$jobs = @()
$seed = $SeedBase
foreach ($fam in $Fams.Split(",")) { foreach ($k in $Ks.Split(",")) { foreach ($d in $Deltas.Split(",")) {
    $seed += 1
    if ($Thetas -ne "") { $th = $Thetas } else { $th = "chain,a+"; if ($fam -eq "random") { $th = "chain,0.02,a+" } }
    $dd = $d.Replace("/","_")
    $outf = ("out/camp" + $Tag + "_" + $fam + "_k" + $k + "_d" + $dd + ".jsonl")
    $jobs += ,@($fam, $k, $seed, $d, $th, $outf)
} } }
Add-Content $log ("{0} campaign {1}: {2} jobs" -f (Get-Date -Format s), $Tag, $jobs.Count)
$running = @()
foreach ($j in $jobs) {
    while (($running | Where-Object { -not $_.P.HasExited }).Count -ge $MaxPar) { Start-Sleep -Seconds 2
        foreach ($r in $running) { if (-not $r.P.HasExited -and ((Get-Date) - $r.T).TotalSeconds -gt 570) {
            Stop-Process -Id $r.P.Id -Force; Add-Content $log ("{0} KILLED pid {1} {2}" -f (Get-Date -Format s), $r.P.Id, $r.Name) } } }
    $ok = $false
    for ($w = 0; $w -lt 15; $w++) {
        $free = [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1024)
        if ($free -gt 800) { $ok = $true; break }
        Add-Content $log ("{0} wait: free {1} MB" -f (Get-Date -Format s), $free); Start-Sleep -Seconds 20
    }
    if (-not $ok) { Add-Content $log ("{0} SKIP (low memory) {1}" -f (Get-Date -Format s), $j[5]); continue }
    $argl = "tbx_run.py {0} {1} {2} {3} {4} {5} {6} {7} {8}" -f $j[5], $j[0], $j[1], $j[2], $NPack, $j[3], $j[4], $NY, $TLimit
    $err = $j[5] + ".err"; $out = $j[5] + ".out"
    $p = Start-Process -FilePath $py -ArgumentList $argl -NoNewWindow -PassThru -RedirectStandardOutput $out -RedirectStandardError $err
    Add-Content $log ("{0} START pid {1} free {2} MB : {3}" -f (Get-Date -Format s), $p.Id, $free, $argl)
    $running += [pscustomobject]@{ P = $p; T = (Get-Date); Name = $j[5] }
}
while (($running | Where-Object { -not $_.P.HasExited }).Count -gt 0) { Start-Sleep -Seconds 2
    foreach ($r in $running) { if (-not $r.P.HasExited -and ((Get-Date) - $r.T).TotalSeconds -gt 570) {
        Stop-Process -Id $r.P.Id -Force; Add-Content $log ("{0} KILLED pid {1} {2}" -f (Get-Date -Format s), $r.P.Id, $r.Name) } } }
foreach ($r in $running) { Add-Content $log ("{0} END pid {1} exit {2} : {3}" -f (Get-Date -Format s), $r.P.Id, $r.P.ExitCode, (Get-Content ($r.Name + ".out") -ErrorAction SilentlyContinue)) }
Add-Content $log ("{0} campaign {1} finished" -f (Get-Date -Format s), $Tag)

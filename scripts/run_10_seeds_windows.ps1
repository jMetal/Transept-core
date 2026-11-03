param(
    [Parameter(Mandatory = $true)]
    [string]$JarPath,

    [Parameter(Mandatory = $true)]
    [string]$DataDirectory,

    [Parameter(Mandatory = $true)]
    [string]$ProblemName,

    [Parameter(Mandatory = $true)]
    [string]$OutputRoot,

    # Comma-separated tools whose precomputed alignments seed the population,
    # e.g. "clustalw,kalign,mafft,muscle5,tcoffee,tmaligner"
    [Parameter(Mandatory = $true)]
    [string]$SeedMethods,

    [int]$MaxEvaluations = 2500,
    [int]$PopulationSize = 100,
    [int]$NumberOfCores = 4,
    [int]$ObserverFrequency = 100
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if (-not (Test-Path -LiteralPath $JarPath -PathType Leaf)) {
    throw "JAR file not found: $JarPath"
}

if (-not (Test-Path -LiteralPath $DataDirectory -PathType Container)) {
    throw "Dataset directory not found: $DataDirectory"
}

if ($MaxEvaluations -le 0) {
    throw "MaxEvaluations must be greater than zero."
}

if ($PopulationSize -le 0) {
    throw "PopulationSize must be greater than zero."
}

if ($NumberOfCores -le 0) {
    throw "NumberOfCores must be greater than zero."
}

if ($ObserverFrequency -le 0) {
    throw "ObserverFrequency must be greater than zero."
}

$Seeds = @(0, 21, 42, 63, 84, 105, 126, 147, 168, 189)

Write-Host "Dataset       : $DataDirectory"
Write-Host "Problem       : $ProblemName"
Write-Host "Output root   : $OutputRoot"
Write-Host "Evaluations   : $MaxEvaluations"
Write-Host "Population    : $PopulationSize"
Write-Host "Cores         : $NumberOfCores"
Write-Host "Seed methods  : $SeedMethods"
Write-Host "Seeds         : $($Seeds -join ', ')"
Write-Host ""

for ($Index = 0; $Index -lt $Seeds.Count; $Index++) {
    $Seed = [long]$Seeds[$Index]
    $RunId = "run_{0:D2}" -f ($Index + 1)

    Write-Host "============================================================"
    Write-Host "Starting $ProblemName / $RunId / seed=$Seed"
    Write-Host "============================================================"

    & java -jar $JarPath run `
        --data-root $DataDirectory `
        --dataset $ProblemName `
        --seeds $SeedMethods `
        --output-root $OutputRoot `
        --run-id $RunId `
        --seed $Seed `
        --max-evaluations $MaxEvaluations `
        --population-size $PopulationSize `
        --threads $NumberOfCores `
        --progress-interval $ObserverFrequency

    if ($LASTEXITCODE -ne 0) {
        throw "Execution failed for $RunId with seed $Seed. Exit code: $LASTEXITCODE"
    }

    $RunDirectory = Join-Path $OutputRoot "ejecuciones\$ProblemName\$RunId"
    $RequiredFiles = @(
        "FUN.csv",
        "progress.jsonl",
        "metadata.json",
        "runtime.csv",
        "seed.txt"
    )

    foreach ($RequiredFile in $RequiredFiles) {
        $RequiredPath = Join-Path $RunDirectory $RequiredFile
        if (-not (Test-Path -LiteralPath $RequiredPath -PathType Leaf)) {
            throw "Expected output file was not generated: $RequiredPath"
        }
    }

    $MsaFiles = @(Get-ChildItem -LiteralPath $RunDirectory -Filter "MSASol*.fasta" -File)
    if ($MsaFiles.Count -eq 0) {
        throw "No MSASol*.fasta files were generated in: $RunDirectory"
    }

    Write-Host "Completed $RunId with seed $Seed"
    Write-Host "Results: $RunDirectory"
    Write-Host ""
}

Write-Host "All ten reproducible executions completed successfully."

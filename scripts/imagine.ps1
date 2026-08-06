<#
.SYNOPSIS
  Generate an image via APINebula (or compatible) using grok-imagine-image.

.DESCRIPTION
  APINebula exposes grok-imagine-image on /v1/chat/completions (not /v1/images/*).
  This script posts a prompt, extracts the image URL from the reply, and downloads it.

.PARAMETER Prompt
  Image description (required unless passed as remaining args).

.PARAMETER OutFile
  Output path. Default: .\images\imagine-yyyyMMdd-HHmmss.jpg

.PARAMETER BaseUrl
  API base URL. Default: https://apinebula.com/v1

.PARAMETER Model
  Model id. Default: grok-imagine-image

.PARAMETER ApiKey
  API key. Default: $env:APINEBULA_API_KEY (fallback $env:XAI_API_KEY)

.PARAMETER NoDownload
  Only print the image URL; do not download.

.EXAMPLE
  .\scripts\imagine.ps1 "iPhone 17 silver product photo, white studio background"

.EXAMPLE
  .\scripts\imagine.ps1 -Prompt "sunset over ocean" -OutFile .\sunset.jpg
#>
[CmdletBinding()]
param(
    [Parameter(Position = 0, ValueFromRemainingArguments = $true)]
    [string[]]$PromptParts,

    [string]$Prompt,

    [string]$OutFile,

    [string]$BaseUrl = "https://apinebula.com/v1",

    [string]$Model = "grok-imagine-image",

    [string]$ApiKey,

    [switch]$NoDownload
)

$ErrorActionPreference = "Stop"

function Get-ApiKey {
    param([string]$Explicit)
    if ($Explicit) { return $Explicit }
    if ($env:APINEBULA_API_KEY) { return $env:APINEBULA_API_KEY }
    if ($env:XAI_API_KEY) { return $env:XAI_API_KEY }
    throw "No API key. Set APINEBULA_API_KEY (or pass -ApiKey)."
}

function Get-PromptText {
    param([string]$Prompt, [string[]]$PromptParts)
    if ($Prompt -and $Prompt.Trim()) { return $Prompt.Trim() }
    if ($PromptParts -and $PromptParts.Count -gt 0) {
        return ($PromptParts -join " ").Trim()
    }
    throw "Missing prompt. Example: .\scripts\imagine.ps1 `"a red apple on a table`""
}

function ConvertFrom-JsonSafe {
    param([string]$Text)
    try {
        return $Text | ConvertFrom-Json
    } catch {
        throw "Failed to parse JSON response: $($_.Exception.Message)`n$Text"
    }
}

function Find-ImageUrl {
    param([string]$Text)
    if (-not $Text) { return $null }

    # Markdown: ![alt](https://...)
    $m = [regex]::Match($Text, '!\[[^\]]*\]\((https?://[^)\s]+)\)')
    if ($m.Success) { return $m.Groups[1].Value }

    # Bare https URL ending with common image extensions
    $m = [regex]::Match($Text, '(https?://[^\s`"''<>]+\.(?:jpg|jpeg|png|webp|gif)(?:\?[^\s`"''<>]*)?)', 'IgnoreCase')
    if ($m.Success) { return $m.Groups[1].Value }

    # Any https URL on known image hosts
    $m = [regex]::Match($Text, '(https?://(?:pubimage|cdn|img)[^\s`"''<>]+)', 'IgnoreCase')
    if ($m.Success) { return $m.Groups[1].Value.TrimEnd(').,;') }

    return $null
}

$promptText = Get-PromptText -Prompt $Prompt -PromptParts $PromptParts
$key = Get-ApiKey -Explicit $ApiKey
$BaseUrl = $BaseUrl.TrimEnd("/")
$endpoint = "$BaseUrl/chat/completions"

if (-not $OutFile) {
    $dir = Join-Path (Get-Location) "images"
    if (-not (Test-Path $dir)) {
        New-Item -ItemType Directory -Path $dir | Out-Null
    }
    $stamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $OutFile = Join-Path $dir "imagine-$stamp.jpg"
}

$bodyObj = @{
    model = $Model
    messages = @(
        @{
            role = "user"
            content = $promptText
        }
    )
}
$bodyJson = $bodyObj | ConvertTo-Json -Depth 6 -Compress
$tmpBody = Join-Path $env:TEMP ("grok-imagine-" + [guid]::NewGuid().ToString("N") + ".json")
[System.IO.File]::WriteAllText($tmpBody, $bodyJson)

Write-Host "POST $endpoint"
Write-Host "model: $Model"
Write-Host "prompt: $promptText"

try {
    $raw = & curl.exe -sS --fail-with-body `
        -X POST $endpoint `
        -H "Authorization: Bearer $key" `
        -H "Content-Type: application/json" `
        --data-binary "@$tmpBody"
    if ($LASTEXITCODE -ne 0) {
        throw "HTTP request failed (exit $LASTEXITCODE): $raw"
    }
} finally {
    Remove-Item -LiteralPath $tmpBody -ErrorAction SilentlyContinue
}

$resp = ConvertFrom-JsonSafe -Text $raw
$content = $null
if ($resp.choices -and $resp.choices[0].message.content) {
    $content = [string]$resp.choices[0].message.content
} elseif ($resp.output_text) {
    $content = [string]$resp.output_text
} else {
    $content = $raw
}

Write-Host ""
Write-Host "assistant:"
Write-Host $content

$url = Find-ImageUrl -Text $content
if (-not $url) {
    throw "No image URL found in model response. Full JSON:`n$raw"
}

Write-Host ""
Write-Host "image url: $url"

if ($NoDownload) {
    Write-Output $url
    exit 0
}

$outDir = Split-Path -Parent $OutFile
if ($outDir -and -not (Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}

& curl.exe -sS -L --fail -o $OutFile $url
if ($LASTEXITCODE -ne 0) {
    throw "Failed to download image from $url"
}

$item = Get-Item -LiteralPath $OutFile
Write-Host ""
Write-Host ("saved: {0} ({1:N0} bytes)" -f $item.FullName, $item.Length)
Write-Output $item.FullName

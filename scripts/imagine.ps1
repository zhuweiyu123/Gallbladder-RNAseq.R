# ========================================================================
# 【中文阅读指南】通过 HTTP 接口生成并下载图片的辅助脚本
# 输入：Prompt 或位置参数中的图片描述；ApiKey 参数或环境变量提供认证信息。
# 流程：整理参数 → 构造 JSON → curl.exe 发送请求 → 解析返回文字中的图片 URL → 下载到 OutFile。
# 输出：默认当前目录 images 下的时间戳 JPG 文件；-NoDownload 只输出链接。
# 这是仓库中的独立图片辅助工具，不参与单细胞分析。BaseUrl 和 Model 控制请求目标。
# 本文件注释解释现有实现；实际接口是否仍支持该模型应在需要调用时另行确认。
# PowerShell 入门：$ 开头是变量；@{} 是键值表；| 把结果传给下一条命令；# 后为注释。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
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
# 【参数入口】运行时可传 -Prompt、-OutFile 等；[switch] 表示只需写开关名，不必再跟 True。
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

# 【函数：Get-ApiKey】按优先级取认证信息：显式 ApiKey → APINEBULA_API_KEY → XAI_API_KEY。
# 都没有时停止；不在脚本中写入实际密钥。
function Get-ApiKey {
    param([string]$Explicit)
    if ($Explicit) { return $Explicit }
    if ($env:APINEBULA_API_KEY) { return $env:APINEBULA_API_KEY }
    if ($env:XAI_API_KEY) { return $env:XAI_API_KEY }
    throw "No API key. Set APINEBULA_API_KEY (or pass -ApiKey)."
}

# 【函数：Get-PromptText】优先使用命名参数 Prompt，否则把位置参数 PromptParts 用空格拼起来。
# 去掉首尾空白；没有有效描述时抛出错误。
function Get-PromptText {
    param([string]$Prompt, [string[]]$PromptParts)
    if ($Prompt -and $Prompt.Trim()) { return $Prompt.Trim() }
    if ($PromptParts -and $PromptParts.Count -gt 0) {
        return ($PromptParts -join " ").Trim()
    }
    throw "Missing prompt. Example: .\scripts\imagine.ps1 `"a red apple on a table`""
}

# 【函数：ConvertFrom-JsonSafe】把接口返回的 JSON 字符串转为 PowerShell 对象，解析失败时提供错误上下文。
function ConvertFrom-JsonSafe {
    param([string]$Text)
    try {
        return $Text | ConvertFrom-Json
    } catch {
        throw "Failed to parse JSON response: $($_.Exception.Message)`n$Text"
    }
}

# 【函数：Find-ImageUrl】依次尝试 Markdown 图片、常见图片扩展名 URL、指定图片主机形式。
# 返回第一个匹配链接；都没有时返回 null，由主流程判断为失败。
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

# 【请求体】哈希表描述模型和用户消息；下一步转换成接口可接收的 JSON 文本。
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
# 【临时文件】将 JSON 写入临时文件，通过 --data-binary 发送，减少命令行引号和中文编码干扰。
$tmpBody = Join-Path $env:TEMP ("grok-imagine-" + [guid]::NewGuid().ToString("N") + ".json")
[System.IO.File]::WriteAllText($tmpBody, $bodyJson)

Write-Host "POST $endpoint"
Write-Host "model: $Model"
Write-Host "prompt: $promptText"

try {
    # 【发送请求】& 调用外部程序，反引号续行；Authorization 携带令牌，返回文本存入 raw。
    $raw = & curl.exe -sS --fail-with-body `
        -X POST $endpoint `
        -H "Authorization: Bearer $key" `
        -H "Content-Type: application/json" `
        --data-binary "@$tmpBody"
    if ($LASTEXITCODE -ne 0) {
        throw "HTTP request failed (exit $LASTEXITCODE): $raw"
    }
# 【清理】请求成功或失败都会进入 finally，删除本次创建的请求体临时文件。
} finally {
    Remove-Item -LiteralPath $tmpBody -ErrorAction SilentlyContinue
}

# 【响应解析】优先读 choices[0].message.content，再尝试 output_text；最后保留原始响应用于排查。
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

# 【仅返回链接】启用 -NoDownload 时此处正常退出，后续下载代码不再执行。
if ($NoDownload) {
    Write-Output $url
    exit 0
}

$outDir = Split-Path -Parent $OutFile
if ($outDir -and -not (Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}

# 【下载图片】-L 跟随重定向，-o 写入 OutFile；退出码非零时报告下载失败。
& curl.exe -sS -L --fail -o $OutFile $url
if ($LASTEXITCODE -ne 0) {
    throw "Failed to download image from $url"
}

$item = Get-Item -LiteralPath $OutFile
Write-Host ""
Write-Host ("saved: {0} ({1:N0} bytes)" -f $item.FullName, $item.Length)
Write-Output $item.FullName

# Text Encoding Rules

Chinese text, pasted notes, Markdown, CSV/TSV metadata, R scripts, and Python templates must be treated as UTF-8.

## Required Defaults

- Read text files with UTF-8 or UTF-8 with BOM support.
- In Python, prefer `encoding="utf-8-sig"` for user-supplied CSV/TSV/text and `encoding="utf-8"` for generated repo files.
- In R, prefer `fileEncoding = "UTF-8-BOM"` for user-supplied CSV files and `encoding = "UTF-8"` for `readLines()`.
- Do not rely on Windows ANSI, locale defaults, or PowerShell default decoding for Chinese text.
- If a file still looks garbled, inspect encoding/provenance first; do not reinterpret the biological content until the text is readable.

## PowerShell Reading

When manually reading Chinese attachments or pasted text in this workspace, use an explicit UTF-8 command such as:

```powershell
Get-Content -LiteralPath <path> -Encoding UTF8
```

For Python-based inspection:

```python
Path(path).read_text(encoding="utf-8-sig")
```

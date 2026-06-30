# PATH Setup

## Required PATH Entries

Add these directories to user PATH when they exist. R must be 4.3 or newer; Python modules use Python 3.10 or newer.

```text
<R_HOME>\bin
E:\rtools44\usr\bin
E:\rtools44\x86_64-w64-mingw32.static.posix\bin
C:\rtools44\usr\bin
C:\rtools44\x86_64-w64-mingw32.static.posix\bin
E:\rtools43\usr\bin
E:\rtools43\x86_64-w64-mingw32.static.posix\bin
C:\rtools43\usr\bin
C:\rtools43\x86_64-w64-mingw32.static.posix\bin
C:\Program Files\Git\usr\bin
C:\ProgramData\miniconda3
C:\ProgramData\miniconda3\Scripts
C:\ProgramData\miniconda3\condabin
```

Also add JAGS if installed, commonly:

```text
C:\Program Files\JAGS\JAGS-4.3.1\x64\bin
C:\Program Files\JAGS\JAGS-4.3.1\bin
```

## Verification Commands

```powershell
Rscript --version
where.exe Rscript
where.exe make
where.exe gcc
where.exe gzip
conda --version
conda run -n seuratv5-course-py python -c "import scanpy, celltypist; print('ok')"
```

## Notes

- Put the selected R 4.3 or newer installation before older Conda R entries so `Rscript` resolves to the correct runtime.
- Rtools must match the installed R minor version, for example Rtools43 for R 4.3 and Rtools44 for R 4.4.
- If package compilation fails, check that `make`, `gcc`, and `g++` resolve from the matching Rtools installation.
- If JAGS-linked packages fail, install JAGS and restart the shell after PATH changes.

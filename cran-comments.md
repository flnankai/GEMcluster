## R CMD check results

This submission includes the arXiv reference for the companion paper:
<https://arxiv.org/abs/2605.08995>.

Local Windows check was run with:

```sh
R CMD check --no-manual GEMcluster_0.1.0.tar.gz
```

Result:

```text
0 errors | 0 warnings | 1 note
```

The remaining NOTE is Windows-specific:

```text
Note: information on .o files for x64 is not available
```

This is platform-related and does not indicate an issue in the package code.
The local check machine had restricted network access, so repository-index
lookup messages may also appear in the console output; they are not package
check failures.

## Downstream dependencies

This is the first CRAN submission. There are no downstream dependencies.

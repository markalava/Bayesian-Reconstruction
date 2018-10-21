# Contraceptive Use All Women package

R package that provides contraceptive use estimates for all women.

Currently exports only the functions required to generate one-country run tables (for integration with FPET). The code to create global runs and produce charts, etc. is included, but not specifically exported.

To install the package, go the the package's parent folder and use `devtools::install`.

```
setwd("..")
devtools::install("FPET-all-women")
```

It can also be installed from GitHub. This requires access to the private GitHub repository and a personal access token (PAT).

```
devtools::install_github("AvenirHealth/FPET-all-women", auth_token = devtools::github_pat("<PAT>"))
```

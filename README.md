# ms-scverse-bioconductor

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/lucas-diedrich/ms-scverse-bioconductor/test.yaml?branch=main
[badge-docs]: https://app.readthedocs.org/projects/ms-scverse-bioconductor/badge/

Joint MS-proteomics analysis of Python and R packages

> [!NOTE]
> Work in progress


This repository documents approaches on how to integrate R and Python analyses

## For users
This GitHub repository contains tutorials on how to perform cross-ecosystem MS-proteomics data analyses in Python and R, specifically with alphapepttools/mulink (in Python) and scp/QFeatures (in R). Please refer to the [documentation][], in particular, the [tutorials][].

## For developers

### Installation

We provide some utility functions to generate test data for developers. You can install this repository as a regular python package with command line utilities.

```
conda create -n msproteomics python=3.13 -y
pip install git+https://github.com/lucas-diedrich/ms-scverse-bioconductor.git
```

### Usage
Use the command line utilities to generate test data etc.

```bash
ms-scverse-bioconductor --help
```

### Docker
If you do not want to install a python environment, you can also use the docker image that mirrors the CLI interface:

```bash
docker run --rm -it ldiedrich/ms-scverse-bioconductor --help
```

## Release notes

See the [changelog][].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][].
If you found a bug, please use the [issue tracker][].


[uv]: https://github.com/astral-sh/uv
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/lucas-diedrich/ms-scverse-bioconductor/issues
[tests]: https://github.com/lucas-diedrich/ms-scverse-bioconductor/actions/workflows/test.yaml
[documentation]: https://github.io/lucas-diedrich/ms-scverse-bioconductor
[changelog]: https://github.com/lucas-diedrich/ms-scverse-bioconductor/releases
[api documentation]: https://ms-scverse-bioconductor.readthedocs.io/page/api.html
[tutorials]: https://ms-scverse-bioconductor.readthedocs.io/page/api.html

[venv]: https://docs.python.org/3/tutorial/venv.html

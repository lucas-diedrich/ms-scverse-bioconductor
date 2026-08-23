# Image that exposes the `ms-scverse-bioconductor` command line interface.
#
# Build (GIT_REF is any branch, tag or commit, GIT_REPOSITORY any fork of the repository):
#   docker build -t ms-scverse-bioconductor .
#   docker build --build-arg GIT_REF=v0.0.1 -t ms-scverse-bioconductor .
#
# Run: every argument is forwarded to the CLI, files are written to /data
#   docker run --rm ms-scverse-bioconductor --help
#   docker run --rm -v "$PWD:/data" ms-scverse-bioconductor export --n-mod 3 --n-obs 50 --name simulation
FROM python:3.14-slim

ARG GIT_REPOSITORY=lucas-diedrich/ms-scverse-bioconductor
ARG GIT_REF=main

# Installing from the source tarball instead of `git+https://…` keeps git out of the image.
RUN pip install --no-cache-dir \
  "ms-scverse-bioconductor @ https://github.com/${GIT_REPOSITORY}/archive/${GIT_REF}.tar.gz"

# The export subcommand writes relative to the working directory, so this is the mount point for results.
WORKDIR /data

ENTRYPOINT [ "ms-scverse-bioconductor" ]
CMD [ "--help" ]

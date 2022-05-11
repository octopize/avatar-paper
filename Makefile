SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

install:  ## Install everything
	poetry install
	Rscript -e 'install.packages(c("IRkernel", "librarian", "styler"), repos="https://cloud.r-project.org")'
	RScript -e 'IRkernel::installspec(user = TRUE)'
.PHONY: install

lint-fix:  ## Format code
	Rscript -e 'library("styler"); style_dir(".")'
.PHONY: lint-fix

notebook:  ## Start jupyter notebook
	poetry run jupyter notebook
.PHONY: notebook


.DEFAULT_GOAL := help
help: Makefile
	@awk 'BEGIN {FS = ":.*##"; printf "Usage: make \033[36m<target>\033[0m\n"} /^[\/\.a-zA-Z1-9_-]+:.*?##/ { printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

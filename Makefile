.PHONY: clean lint format test coverage build develop
.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: ## remove build artifacts
	@rm -rf build/
	@rm -rf dist/

lint:  ## Check for python formatting issues via black & flake8
	@ruff format . --check && ruff check .

format:  ## Modify python code using black & show flake8 issues
	@ruff format . && ruff check . --fix --show-fixes

test: ## Run unit tests
	@py.test

coverage: ## Generate coverage report
	@coverage run -m pytest
	@coverage html

build: ## Rebuild in development environment
	@python setup.py develop
	@ls -lh src/pybmds

dist: clean ## Build wheel package for distribution
	@python setup.py bdist_wheel
	@ls -l dist

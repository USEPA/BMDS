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
	@rm -rf build/ dist/

lint:  ## Check for python formatting issues
	@ruff format . --check && ruff check .

format:  ## Modify python code where possible
	@ruff format . && ruff check . --fix --show-fixes

test: ## Run unit tests
	@py.test

coverage: ## Generate coverage report
	@coverage run -m pytest
	@coverage html

build: ## Rebuild in development environment
	@python setup.py develop
	@stubgen -p pybmds.bmdscore -o src/
	@ruff format src/pybmds/bmdscore.pyi
	@ls -lh src/pybmds/

dist: ## Build wheel package for distribution
	@python setup.py bdist_wheel
	@ls -lh dist

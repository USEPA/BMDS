.PHONY: clean lint format test coverage build
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
	@black src/pybmds tests --check && flake8 src/pybmds tests

format:  ## Modify python code using black & show flake8 issues
	@black src/pybmds tests && isort src/pybmds tests && flake8 src/pybmds tests

test: ## Run unit tests
	@py.test

coverage: ## Generate coverage report
	@coverage run -m pytest
	@coverage html
	@open htmlcov/index.html

build: clean ## Build wheel package
	@python setup.py bdist_wheel
	@ls -l dist

.PHONY: clean lint format test build
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
	@black . --check && flake8 .

format:  ## Modify python code using black & show flake8 issues
	@black . && isort . && flake8 .

test:
	@py.test

build: clean ## build wheel package
	@python setup.py bdist_wheel
	@ls -l dist

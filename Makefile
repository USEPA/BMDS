.PHONY: clean lint format format-cpp test coverage build dist docs docs-docx docs-serve docs-clean loc
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

format-cpp: ## Format C++ code
	@find src -name "*.cpp" -o -name "*.h" | xargs clang-format -i

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

docs: ## Build documentation {html}
	sphinx-build -W -b html docs/source docs/build/html
	@echo "HTML: \"docs/build/html/index.html\""

docs-docx: ## Build documentation {docx}
	sphinx-build -W -b singlehtml docs/source docs/build/singlehtml
	cd docs/build/singlehtml && pandoc -s index.html -o ../pybmds.docx
	@echo "docx: \"docs/build/pybmds.docx\""

docs-serve: ## Realtime documentation preview
	sphinx-autobuild -b html docs/source docs/build/html --port 5800

docs-clean: ## Clean documentation
	@$(MAKE) -C docs clean

loc: ## Generate lines of code report
	@cloc \
		--exclude-dir=build,dist,venv \
		--exclude-ext=json,yaml,svg,toml,ini \
		--vcs=git \
		--counted loc-files.txt \
		--md \
		.

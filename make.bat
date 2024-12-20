@ECHO off

if "%~1" == "" goto :help
if /I %1 == help goto :help
if /I %1 == docs goto :docs
if /I %1 == docs-serve goto :docs-serve
if /I %1 == lint goto :lint
if /I %1 == format goto :format
if /I %1 == test goto :test
if /I %1 == coverage goto :coverage
if /I %1 == build goto :build
if /I %1 == dist goto :dist
goto :help

:help
echo.Please use `make ^<target^>` where ^<target^> is one of
echo.  docs         Build documentation {html}
echo.  docs-serve   Realtime documentation preview
echo.  lint         perform both lint-py and lint-js
echo.  format       perform both format-py and lint-js
echo.  test         run python tests
echo.  coverage     generate coverage report
echo.  build        rebuild in development environment
echo.  dist         build wheel package for distribution
goto :eof

:docs
rmdir /s /q docs\build
sphinx-build -W -b html docs/source docs/build/html
goto :eof

:docs-serve
rmdir /s /q docs\build
sphinx-autobuild -b html docs/source docs/build/html --port 5800
goto :eof

:lint
ruff format . --check && ruff check .
goto :eof

:format
ruff format . && ruff check . --fix --show-fixes
goto :eof

:test
py.test
goto :eof

:coverage
coverage run -m pytest
coverage html
goto :eof

:build
python setup.py develop
stubgen -p pybmds.bmdscore -o src
ruff format src\pybmds\bmdscore.pyi
goto :eof

:dist
rmdir /s /q .\build
rmdir /s /q .\dist
python setup.py bdist_wheel
dir .\dist
goto :eof

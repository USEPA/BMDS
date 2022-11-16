@ECHO off

if "%~1" == "" goto :help
if /I %1 == help goto :help
if /I %1 == lint goto :lint
if /I %1 == format goto :format
if /I %1 == test goto :test
if /I %1 == build goto :build
goto :help

:help
echo.Please use `make ^<target^>` where ^<target^> is one of
echo.  test         run python tests
echo.  test-refresh removes mock requests and runs python tests
echo.  lint         perform both lint-py and lint-js
echo.  format       perform both format-py and lint-js
echo.  build        build a python wheel for distribution
goto :eof

:lint
black . --check && flake8 .
goto :eof

:format
black . && isort . && flake8 .
goto :eof

:test
py.test
goto :eof

:build
rmdir /s /q .\build
rmdir /s /q .\dist
python setup.py bdist_wheel
dir .\dist
goto :eof

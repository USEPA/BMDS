set -x

# pip install -e .
for run in {1..10}; do
  py.test tests/test_pybmds/tests/test_sessions.py -k test_nested_dichotomous -s
done

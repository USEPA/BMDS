from pybmds.models import nested_dichotomous


def test_execution(nd_dataset4):
    # add seed for reproducibility
    analysis = nested_dichotomous.NestedLogistic(nd_dataset4, settings=dict(bootstrap_seed=1))
    analysis.execute()
    text = analysis.text()
    assert len(text) > 0

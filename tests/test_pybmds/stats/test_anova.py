from pybmds.stats.anova import AnovaTests


def test_output_when_anova_cannot_be_calculated():
    assert AnovaTests.output_3tests(None) == "ANOVA cannot be calculated for this dataset."

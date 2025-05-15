import pytest
from extra_nested_plots import generate_extra_plots

import pybmds


@pytest.mark.mpl_image_compare
def test_nested_plot(nd_dataset4):
    session = pybmds.Session(dataset=nd_dataset4)
    session.add_default_models(
        settings={"litter_specific_covariate": 1, "intralitter_correlation": 1}
    )
    session.execute()
    return generate_extra_plots(session)

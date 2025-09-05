# import numpy as np
import pytest

from pybmds.models import continuous


class TestBmdModelContinuousIndividual:
    def test_increasing(self, cidataset):
        for Model in [
            continuous.ExponentialM3,
            continuous.ExponentialM5,
            continuous.Power,
            continuous.Hill,
            continuous.Linear,
            continuous.Polynomial,
        ]:
            result = Model(cidataset).execute()
            assert result.has_completed and result.bmd > 0, Model.__name__

    @pytest.mark.mpl_image_compare
    def test_continuous_individual_plot(self, cidataset):
        model = continuous.Power(dataset=cidataset)
        model.execute()
        return model.plot()

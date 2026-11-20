import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from proteomicscopylot import notebook_functions as nb


@pytest.mark.parametrize("plot", [nb.boxplot_int, nb.histo_int])
def test_linear_and_log2_inputs_render_identically(plot):
    meta = pd.DataFrame({"sample": list("abcd"), "condition": ["A", "A", "B", "B"]})
    logged = pd.DataFrame([[0, 4, 2, 1], [-1, 1, 3, 5], [4, 2, 0, -2],
                           [np.nan, 3, 4, 5], [6, 0, 1, np.nan]], columns=list("abcd"))
    linear = np.exp2(logged)
    linear.loc[3, "a"] = 0  # A raw non-detection; linear 1 above is a measured log2 zero.

    def pixels(frame, **kwargs):
        figure = plot(frame, meta, **kwargs)
        try:
            figure.canvas.draw()
            return np.asarray(figure.canvas.buffer_rgba()).copy()
        finally:
            plt.close(figure)

    expected = pixels(logged)
    np.testing.assert_array_equal(pixels(logged, data_scale="log2"), expected)
    np.testing.assert_array_equal(pixels(linear, data_scale="linear"), expected)
    with pytest.raises(ValueError, match="data_scale"):
        plot(linear, meta, data_scale="unknown")

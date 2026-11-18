"""Crowded legends must preserve both readable keys and a useful data area."""
import pytest


@pytest.mark.parametrize('count', [1, 2, 5, 10, 20, 40])
@pytest.mark.parametrize('long_names', [False, True])
def test_dynamic_legend_preserves_keys_styles_and_data(count, long_names):
    import matplotlib.pyplot as plt
    import numpy as np
    from app.services.plot_export import to_png_bytes

    fig, ax = plt.subplots(figsize=(24 / 2.54, 16 / 2.54))
    labels = [f'Condition {i:02}' + (' long treatment name' if long_names else '') for i in range(count)]
    for i, label in enumerate(labels):
        ax.scatter([0, 1], [i, i + 1], label=label)
    old = ax.legend(title='Treatment', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)
    before = [c.get_offsets().copy() for c in ax.collections]
    colors = [h.get_facecolors().copy() for h in old.legend_handles]
    ax.set_xlabel('Abundance rank')
    ax.set_ylabel('Intensity')
    ax.set_title('Synthetic abundance')
    fig.tight_layout(rect=[0, 0, .85, 1])
    assert to_png_bytes(fig, plt).startswith(b'\x89PNG')
    legend = ax.get_legend()
    assert [t.get_text() for t in legend.get_texts()] == labels
    assert [t.get_fontsize() for t in legend.get_texts()] == [10] * count
    assert legend.get_title().get_text() == 'Treatment'
    for actual, expected in zip(ax.collections, before):
        np.testing.assert_array_equal(actual.get_offsets(), expected)
    for actual, expected in zip(legend.legend_handles, colors):
        np.testing.assert_array_equal(actual.get_facecolors(), expected)
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig)
    fig.canvas.draw()
    extent = legend.get_window_extent(fig.canvas.get_renderer())
    assert extent.x0 >= 0 and extent.y0 >= 0
    assert extent.x1 <= fig.bbox.width and extent.y1 <= fig.bbox.height
    if count >= 20:
        assert ax.bbox.height / fig.dpi >= 3.5
        assert ax.bbox.width / fig.dpi >= 6
        assert extent.y1 < ax.bbox.y0
    else:
        assert legend is old


def test_explicit_legend_offset_is_not_reflowed():
    import matplotlib.pyplot as plt
    from app.services.plot_export import _fit_crowded_legend
    from app.services.plot_styling import _style_legend

    fig, ax = plt.subplots()
    for i in range(40):
        ax.plot([0, 1], [0, i], label=str(i))
    legend = ax.legend()
    _style_legend(legend, {'legendXOffset': 12})
    fig.canvas.draw()
    try:
        _fit_crowded_legend(fig, 150)
        assert ax.get_legend() is legend
        assert not hasattr(fig, '_copylot_legend_bottom')
    finally:
        plt.close(fig)

import shiba2corr


def test_version_is_string():
	assert isinstance(shiba2corr.__version__, str)
	assert shiba2corr.__version__.startswith("v")


def test_public_api_exports():
	for name in (
		"EVENT_TYPES",
		"load_union_event_data",
		"compute_dse_overlap_by_direction",
		"compute_weighted_correlation_dpsi_union",
		"save_weighted_correlation_results",
		"save_scatter_plot_data",
		"plot_weighted_correlation_scatter",
		"plot_dse_overlap_venns",
	):
		assert hasattr(shiba2corr, name), f"missing public symbol: {name}"

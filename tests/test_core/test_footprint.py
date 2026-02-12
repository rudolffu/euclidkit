from euclidkit.core.footprint import get_moc_path, list_available_surveys


def test_list_available_surveys_dr1_includes_cgv_input():
    surveys = list_available_surveys("DR1")
    assert "WIDE" in surveys
    assert "DEEP" in surveys
    assert "BOTH" in surveys
    assert "CGV_INPUT" in surveys


def test_get_moc_path_resolves_cgv_alias():
    moc_path = get_moc_path(survey="cgv", data_release="dr1")
    assert moc_path.name == "cgv_map_dr1input_o13_moc.fits"
    assert moc_path.exists()


def test_get_moc_path_resolves_union_alias():
    moc_path = get_moc_path(survey="union", data_release="DR1")
    assert moc_path.name == "dr1_mer_wide_deep_union_o13_moc.fits"
    assert moc_path.exists()


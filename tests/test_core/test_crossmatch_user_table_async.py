"""Tests for remote user-table crossmatch async execution behavior."""

import json
from unittest.mock import Mock

import numpy as np
from astropy.table import Table

from euclidkit.core.data_access import EuclidArchive


def test_crossmatch_user_table_uses_async_by_default():
    archive = EuclidArchive(environment="IDR")
    archive._logged_in = True
    archive.euclid = Mock()

    # Resolve table name + columns without remote calls.
    archive._resolve_user_table_reference = Mock(return_value="user_test.sample_table")
    archive._get_remote_table_columns = Mock(
        return_value=["sample_table_oid", "source_id", "ra_obj", "dec_obj"]
    )
    preflight_job = Mock()
    preflight_job.get_results.return_value = Table(
        {"min_oid": [1], "max_oid": [100], "n_rows": [100]}
    )
    archive.euclid.launch_job.return_value = preflight_job

    survey_job = Mock()
    survey_job.jobid = "survey-job"
    survey_job.get_results.return_value = Table({"source_id": [1], "object_id": [1], "origin": ["survey"]})
    mode_job = Mock()
    mode_job.jobid = "mode-job"
    mode_job.get_results.return_value = Table({"source_id": [1], "object_id": [1], "origin": ["mode"]})
    archive.euclid.launch_job_async.side_effect = [survey_job, mode_job]

    result = archive.crossmatch_user_table(
        user_table_name="sample_table",
        use_object_id=True,
        full_async=False,
        idr_field="WIDE",
    )

    assert archive.euclid.launch_job_async.call_count == 2
    archive.euclid.launch_job.assert_called_once()
    assert isinstance(result, Table)
    assert len(result) == 1
    assert result["origin"][0] == "survey"

    queries = [call.args[0] for call in archive.euclid.launch_job_async.call_args_list]
    assert "JOIN catalogue.mer_catalogue_wide_survey AS m" in queries[0]
    assert "JOIN catalogue.mer_catalogue_wide_mode AS m" in queries[1]
    assert "NOT EXISTS" in queries[1]
    assert "catalogue.mer_catalogue_wide_survey AS s" in queries[1]
    for field in [
        "det_quality_flag",
        "parent_id",
        "spurious_flag",
        "vis_det",
        "flag_vis",
        "flag_y",
        "flag_j",
        "flag_h",
    ]:
        assert f"m.{field} AS {field}" in queries[0]


def test_crossmatch_user_table_deep_both_queries_both_partitions():
    archive = EuclidArchive(environment="IDR")
    archive._logged_in = True
    archive.euclid = Mock()

    archive._resolve_user_table_reference = Mock(return_value="user_test.sample_table")
    archive._get_remote_table_columns = Mock(
        return_value=["sample_table_oid", "source_id", "ra_obj", "dec_obj"]
    )
    preflight_job = Mock()
    preflight_job.get_results.return_value = Table(
        {"min_oid": [1], "max_oid": [100], "n_rows": [100]}
    )
    archive.euclid.launch_job.return_value = preflight_job

    survey_job = Mock()
    survey_job.jobid = "deep-survey-job"
    survey_job.get_results.return_value = Table(
        {"source_id": [1], "object_id": [1], "origin": ["deep_survey"]}
    )
    mode_job = Mock()
    mode_job.jobid = "deep-mode-job"
    mode_job.get_results.return_value = Table(
        {"source_id": [2], "object_id": [2], "origin": ["deep_mode"]}
    )
    archive.euclid.launch_job_async.side_effect = [survey_job, mode_job]

    meta = archive.crossmatch_user_table(
        user_table_name="sample_table",
        use_object_id=True,
        full_async=True,
        idr_field="DEEP",
        idr_deep_partition="both",
    )

    assert archive.euclid.launch_job_async.call_count == 2
    assert meta["results_downloaded"] is True
    assert meta["result_row_count"] == 2
    assert meta["idr_field"] == "DEEP"
    assert meta["idr_deep_partition"] == "both"
    assert meta["mer_tables"] == [
        "catalogue.mer_catalogue_deep_survey",
        "catalogue.mer_catalogue_deep_mode",
    ]

    queries = [call.args[0] for call in archive.euclid.launch_job_async.call_args_list]
    assert "JOIN catalogue.mer_catalogue_deep_survey AS m" in queries[0]
    assert "JOIN catalogue.mer_catalogue_deep_mode AS m" in queries[1]
    assert "NOT EXISTS" not in queries[0]
    assert "NOT EXISTS" not in queries[1]


def test_crossmatch_user_table_drop_empty_columns_single_async(tmp_path):
    archive = EuclidArchive(environment="IDR")
    archive._logged_in = True
    archive.euclid = Mock()

    archive._resolve_user_table_reference = Mock(return_value="user_test.sample_table")
    archive._get_remote_table_columns = Mock(
        return_value=["sample_table_oid", "source_id", "ra_obj", "dec_obj"]
    )
    preflight_job = Mock()
    preflight_job.get_results.return_value = Table(
        {"min_oid": [1], "max_oid": [100], "n_rows": [100]}
    )
    archive.euclid.launch_job.return_value = preflight_job

    job = Mock()
    job.jobid = "remote-job"
    job.get_results.return_value = Table({
        "source_id": [1, 2],
        "object_id": [101, 102],
        "all_nan": [np.nan, np.nan],
        "zero": [0, 0],
    })
    archive.euclid.launch_job_async.return_value = job
    output_path = tmp_path / "remote_pruned.fits"

    meta = archive.crossmatch_user_table(
        user_table_name="sample_table",
        use_object_id=True,
        full_async=True,
        output_file=output_path,
        drop_empty_columns=True,
    )

    saved = Table.read(output_path, format="fits")
    assert "all_nan" not in saved.colnames
    assert "zero" in saved.colnames
    assert meta["dropped_empty_columns"] == ["all_nan"]
    assert meta["dropped_empty_column_count"] == 1


def test_crossmatch_user_table_drop_empty_columns_chunked_final_only(tmp_path):
    archive = EuclidArchive(environment="PDR")
    archive._logged_in = True
    archive.euclid = Mock()

    archive._resolve_user_table_reference = Mock(return_value="user_test.sample_table")
    archive._get_remote_table_columns = Mock(
        return_value=["sample_table_oid", "source_id", "ra_obj", "dec_obj"]
    )
    preflight_job = Mock()
    preflight_job.get_results.return_value = Table(
        {"min_oid": [1], "max_oid": [3], "n_rows": [2_000_000]}
    )
    archive.euclid.launch_job.return_value = preflight_job

    job = Mock()
    job.jobid = "remote-chunk-job"
    job.get_results.return_value = Table({
        "source_id": [1],
        "object_id": [101],
        "all_nan": [np.nan],
        "zero": [0],
    })
    archive.euclid.launch_job_async.return_value = job
    output_path = tmp_path / "remote_chunk_pruned.fits"

    meta = archive.crossmatch_user_table(
        user_table_name="sample_table",
        use_object_id=True,
        full_async=True,
        output_file=output_path,
        drop_empty_columns=True,
    )

    saved = Table.read(output_path, format="fits")
    part_files = sorted(tmp_path.glob("remote_chunk_pruned_part_*.fits"))
    manifest = json.loads((tmp_path / "remote_chunk_pruned.fits.manifest.json").read_text())

    assert "all_nan" not in saved.colnames
    assert "zero" in saved.colnames
    assert "all_nan" in Table.read(part_files[0], format="fits").colnames
    assert meta["dropped_empty_columns"] == ["all_nan"]
    assert manifest["dropped_empty_columns"] == ["all_nan"]

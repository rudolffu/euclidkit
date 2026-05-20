"""Tests for remote user-table crossmatch async execution behavior."""

from unittest.mock import Mock

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

"""Tests for secure login flows and mag conversion in crossmatch."""

import os
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy.table import Table

from euclidkit.core.data_access import EuclidArchive


class TestSecureLogin:
    def setup_method(self):
        self.archive = EuclidArchive(environment='REG')
        # replace underlying client with a mock
        self.archive.euclid = Mock()

    def test_login_with_credentials_file(self, tmp_path):
        cred = tmp_path / "cred.txt"
        cred.write_text("user\npass\n")
        self.archive.login(credentials_file=str(cred))
        # login should be called with credentials_file argument
        self.archive.euclid.login.assert_called_once_with(credentials_file=str(cred))

    def test_login_with_env_vars(self, monkeypatch):
        monkeypatch.setenv("EUCLID_USER", "env_user")
        monkeypatch.setenv("EUCLID_PASSWORD", "env_pass")
        self.archive.login()
        self.archive.euclid.login.assert_called_once_with(user="env_user", password="env_pass")

    def test_login_with_keyring(self, monkeypatch):
        # Only user in env, password from keyring
        monkeypatch.setenv("EUCLID_USER", "kr_user")
        fake_keyring = Mock()
        fake_keyring.get_password.return_value = "kr_pass"
        with patch("euclidkit.core.data_access.keyring", fake_keyring):
            self.archive.login()  # use_keyring defaults True
        self.archive.euclid.login.assert_called_once_with(user="kr_user", password="kr_pass")

    def test_login_with_prompt(self, monkeypatch):
        # Provide user explicitly; simulate interactive prompt
        fake_getpass = Mock(return_value="typed_pw")
        with patch("euclidkit.core.data_access.getpass", fake_getpass):
            self.archive.login(user="prompt_user", prompt=True)
        self.archive.euclid.login.assert_called_once_with(user="prompt_user", password="typed_pw")


class TestCrossmatchMagnitudes:
    def setup_method(self):
        self.archive = EuclidArchive(environment='REG')
        self.archive.euclid = Mock()

    def test_crossmatch_batch_adds_magnitudes(self, tmp_path):
        # Prepare a small input batch
        batch = Table({
            'ra': [150.0],
            'dec': [2.0],
        })

        # Mock job result returned by archive
        flux_vals = {
            'flux_y_templfit': np.array([10.0]),
            'flux_h_templfit': np.array([20.0]),
            'flux_j_templfit': np.array([30.0]),
            'flux_vis_psf': np.array([40.0]),
        }
        result_table = Table({
            **flux_vals,
            'separation_deg': np.array([1.0 / 3600.0]),
        })

        mock_job = Mock()
        mock_job.get_results.return_value = result_table
        self.archive.euclid.launch_job.return_value = mock_job
        self.archive.euclid.launch_job_async.return_value = mock_job

        # Call the private batch method directly
        out = self.archive._crossmatch_batch(
            batch=batch,
            ra_col='ra',
            dec_col='dec',
            radius=1.0,
            mer_table='catalogue.mer_catalogue'
        )

        # separation converted
        assert 'separation_arcsec' in out.colnames
        assert 'separation_deg' not in out.colnames

        # AB mags added
        for flux_col, mag_col in [
            ('flux_y_templfit', 'mag_y_templfit'),
            ('flux_h_templfit', 'mag_h_templfit'),
            ('flux_j_templfit', 'mag_j_templfit'),
            ('flux_vis_psf', 'mag_vis_psf'),
        ]:
            assert mag_col in out.colnames
            expected = -2.5 * np.log10(result_table[flux_col]) + 23.9
            np.testing.assert_allclose(out[mag_col], expected)


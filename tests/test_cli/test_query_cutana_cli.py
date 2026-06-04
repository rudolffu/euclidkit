"""Tests for query-cutana CLI command."""

import os
import tempfile
from unittest.mock import Mock, patch

import pandas as pd
from astropy.table import Table
from click.testing import CliRunner

from euclidkit.cli.crossmatch_cli import query_cutana


class TestQueryCutanaCLI:
    """Test cases for query-cutana command."""

    def setup_method(self):
        self.runner = CliRunner()
        self.sample_table = Table(
            {
                'object_id': [100001, 100002],
            }
        )

    def _create_temp_sources_file(self):
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
        self.sample_table.write(temp_file.name, format='fits', overwrite=True)
        return temp_file.name

    def test_query_cutana_basic_usage(self):
        """Command should load table, run generator, and report output."""
        sources_file = self._create_temp_sources_file()

        try:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
                mock_archive = Mock()
                mock_generator = Mock()
                mock_generator.archive = mock_archive
                mock_generator.generate_cutana_input.return_value = pd.DataFrame(
                    {'SourceID': ['OBJ_100001', 'OBJ_100002']}
                )

                with patch('euclidkit.core.cutouts.CutoutGenerator', return_value=mock_generator):
                    with patch('euclidkit.utils.io.load_table', return_value=self.sample_table):
                        result = self.runner.invoke(
                            query_cutana,
                            [
                                '--sources',
                                sources_file,
                                '--output',
                                output_file.name,
                                '--instrument',
                                'NISP',
                                '--nisp-filters',
                                'NIR_Y,NIR_H',
                                '--cutout-size',
                                'arcsec',
                                '--cutout-size-value',
                                '15',
                            ],
                        )

                assert result.exit_code == 0
                assert "Cutana query completed: 2 sources with mosaic matches" in result.output
                assert "Cutana input file saved to:" in result.output

                mock_archive.login.assert_called_once()
                mock_archive.logout.assert_called_once()
                call_kwargs = mock_generator.generate_cutana_input.call_args.kwargs
                assert call_kwargs['instrument_name'] == 'NISP'
                assert call_kwargs['nisp_filters'] == ['NIR_Y', 'NIR_H']
                assert call_kwargs['cutout_size'] == 'arcsec'
                assert call_kwargs['cutout_size_value'] == 15.0
                assert call_kwargs['idr_field'] is None
                assert call_kwargs['idr_deep_partition'] == 'survey'
        finally:
            os.unlink(sources_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_cutana_idr_field_is_forwarded(self):
        """IDR field selection should be forwarded when environment is IDR."""
        sources_file = self._create_temp_sources_file()

        try:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
                mock_archive = Mock()
                mock_generator = Mock()
                mock_generator.archive = mock_archive
                mock_generator.generate_cutana_input.return_value = pd.DataFrame(
                    {'SourceID': ['OBJ_100001']}
                )

                with patch('euclidkit.core.cutouts.CutoutGenerator', return_value=mock_generator):
                    with patch('euclidkit.utils.io.load_table', return_value=self.sample_table):
                        result = self.runner.invoke(
                            query_cutana,
                            [
                                '--sources',
                                sources_file,
                                '--output',
                                output_file.name,
                                '--environment',
                                'IDR',
                                '--idr-field',
                                'DEEP',
                            ],
                        )

                assert result.exit_code == 0
                call_kwargs = mock_generator.generate_cutana_input.call_args.kwargs
                assert call_kwargs['idr_field'] == 'DEEP'
                assert call_kwargs['idr_deep_partition'] == 'survey'
        finally:
            os.unlink(sources_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_cutana_idr_deep_partition_is_forwarded(self):
        """IDR DEEP partition selection should be forwarded to Cutana generation."""
        sources_file = self._create_temp_sources_file()

        try:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
                mock_archive = Mock()
                mock_generator = Mock()
                mock_generator.archive = mock_archive
                mock_generator.generate_cutana_input.return_value = pd.DataFrame(
                    {'SourceID': ['OBJ_100001']}
                )

                with patch('euclidkit.core.cutouts.CutoutGenerator', return_value=mock_generator):
                    with patch('euclidkit.utils.io.load_table', return_value=self.sample_table):
                        result = self.runner.invoke(
                            query_cutana,
                            [
                                '--sources',
                                sources_file,
                                '--output',
                                output_file.name,
                                '--environment',
                                'IDR',
                                '--idr-field',
                                'DEEP',
                                '--idr-deep-partition',
                                'both',
                            ],
                        )

                assert result.exit_code == 0
                call_kwargs = mock_generator.generate_cutana_input.call_args.kwargs
                assert call_kwargs['idr_field'] == 'DEEP'
                assert call_kwargs['idr_deep_partition'] == 'both'
        finally:
            os.unlink(sources_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_cutana_rejects_nisp_filters_for_vis(self):
        """NISP filters should fail fast when instrument is VIS."""
        sources_file = self._create_temp_sources_file()

        try:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
                mock_archive = Mock()
                mock_generator = Mock()
                mock_generator.archive = mock_archive

                with patch('euclidkit.core.cutouts.CutoutGenerator', return_value=mock_generator):
                    with patch('euclidkit.utils.io.load_table', return_value=self.sample_table):
                        result = self.runner.invoke(
                            query_cutana,
                            [
                                '--sources',
                                sources_file,
                                '--output',
                                output_file.name,
                                '--instrument',
                                'VIS',
                                '--nisp-filters',
                                'NIR_Y',
                            ],
                        )

                assert result.exit_code == 1
                assert "Error generating Cutana input: --nisp-filters requires --instrument NISP" in result.output
                mock_archive.logout.assert_called_once()
                mock_generator.generate_cutana_input.assert_not_called()
        finally:
            os.unlink(sources_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

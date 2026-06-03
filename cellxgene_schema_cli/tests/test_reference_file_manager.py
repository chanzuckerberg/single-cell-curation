"""Tests for reference_file_manager module."""

from unittest.mock import patch

import pytest
import yaml
from cellxgene_schema.reference_file_manager import ReferenceFileManager


@pytest.fixture
def sample_config(tmp_path):
    """Create a sample reference files config."""
    config = {
        "guidescan_indexes": {
            "human": {
                "organism_id": "NCBITaxon:9606",
                "url": "https://example.com/human_index.tar.gz",
                "description": "Human GuideScan2 index",
            },
            "mouse": {
                "organism_id": "NCBITaxon:10090",
                "url": "https://example.com/mouse_index.tar.gz",
                "description": "Mouse GuideScan2 index",
            },
        },
        "other_files": {
            "plain_file": {
                "url": "https://example.com/data.txt",
                "description": "A plain text file",
            },
            "zip_file": {
                "url": "https://example.com/archive.zip",
                "description": "A zip archive",
            },
        },
    }
    config_path = tmp_path / "reference_files.yml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


@pytest.fixture
def manager(tmp_path, sample_config):
    """Create a ReferenceFileManager with sample config."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    return ReferenceFileManager(str(cache_dir), str(sample_config))


class TestReferenceFileManagerInit:
    """Tests for ReferenceFileManager initialization."""

    def test_init_creates_manager(self, manager):
        """Test that manager is created with config loaded."""
        assert manager.config is not None
        assert "guidescan_indexes" in manager.config

    def test_init_invalid_config_path(self, tmp_path):
        """Test that FileNotFoundError is raised for invalid config."""
        with pytest.raises(FileNotFoundError):
            ReferenceFileManager(str(tmp_path / "cache"), str(tmp_path / "nonexistent.yml"))


class TestGetFileConfig:
    """Tests for _get_file_config method."""

    def test_get_file_config_success(self, manager):
        """Test getting a valid file config."""
        config = manager._get_file_config("guidescan_indexes", "human")
        assert config["organism_id"] == "NCBITaxon:9606"
        assert "url" in config

    def test_get_file_config_invalid_category(self, manager):
        """Test error for invalid category."""
        with pytest.raises(KeyError, match="Category 'invalid_category' not found"):
            manager._get_file_config("invalid_category", "human")

    def test_get_file_config_invalid_key(self, manager):
        """Test error for invalid key."""
        with pytest.raises(KeyError, match="Key 'invalid_key' not found"):
            manager._get_file_config("guidescan_indexes", "invalid_key")


class TestGetKeyByOrganismId:
    """Tests for get_key_by_organism_id method."""

    def test_get_key_found(self, manager):
        """Test finding key by organism ID."""
        key = manager.get_key_by_organism_id("guidescan_indexes", "NCBITaxon:9606")
        assert key == "human"

    def test_get_key_not_found(self, manager):
        """Test None returned for unknown organism ID."""
        key = manager.get_key_by_organism_id("guidescan_indexes", "NCBITaxon:99999")
        assert key is None

    def test_get_key_invalid_category(self, manager):
        """Test None returned for invalid category."""
        key = manager.get_key_by_organism_id("invalid_category", "NCBITaxon:9606")
        assert key is None


class TestFetch:
    """Tests for fetch method."""

    def test_fetch_tar_gz_uses_untar_processor(self, manager):
        """Test that .tar.gz files use Untar processor."""
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            mock_retrieve.return_value = ["/path/to/extracted/file1", "/path/to/extracted/file2"]

            result = manager.fetch("guidescan_indexes", "human")

            # Verify Untar processor was used
            call_kwargs = mock_retrieve.call_args
            processor = call_kwargs.kwargs.get("processor") or call_kwargs[1].get("processor")
            assert processor is not None
            assert type(processor).__name__ == "Untar"

            # Verify result is a list
            assert isinstance(result, list)
            assert len(result) == 2

    def test_fetch_zip_uses_unzip_processor(self, manager):
        """Test that .zip files use Unzip processor."""
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            mock_retrieve.return_value = ["/path/to/extracted/file"]

            manager.fetch("other_files", "zip_file")

            call_kwargs = mock_retrieve.call_args
            processor = call_kwargs.kwargs.get("processor") or call_kwargs[1].get("processor")
            assert processor is not None
            assert type(processor).__name__ == "Unzip"

    def test_fetch_plain_file_no_processor(self, manager):
        """Test that non-archive files don't use a processor."""
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            mock_retrieve.return_value = "/path/to/file.txt"

            manager.fetch("other_files", "plain_file")

            call_kwargs = mock_retrieve.call_args
            processor = call_kwargs.kwargs.get("processor") or call_kwargs[1].get("processor")
            assert processor is None

    def test_fetch_normalizes_string_to_list(self, manager):
        """Test that single string result is normalized to list.

        This tests the bug fix for pooch.retrieve returning str (not list)
        when processor=None, which would cause character-by-character iteration.
        """
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            # Simulate pooch returning a single string (no processor case)
            mock_retrieve.return_value = "/path/to/file.txt"

            result = manager.fetch("other_files", "plain_file")

            # Result should be a list, not a string
            assert isinstance(result, list), "fetch() should always return a list"
            assert result == ["/path/to/file.txt"]

            # Verify iterating doesn't give characters
            for item in result:
                assert item == "/path/to/file.txt"
                assert len(item) > 1  # Not a single character

    def test_fetch_preserves_list_from_processor(self, manager):
        """Test that list result from processor is preserved."""
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            # Simulate pooch returning a list (with processor)
            mock_retrieve.return_value = ["/path/file1.txt", "/path/file2.txt"]

            result = manager.fetch("guidescan_indexes", "human")

            assert isinstance(result, list)
            assert result == ["/path/file1.txt", "/path/file2.txt"]


class TestClearCache:
    """Tests for clear_cache method."""

    def test_clear_cache_all(self, manager, tmp_path):
        """Test clearing all cached files."""
        # Create some cached files
        cache_dir = manager.cache_dir
        (cache_dir / "category1").mkdir()
        (cache_dir / "category1" / "file1.txt").touch()
        (cache_dir / "category2").mkdir()
        (cache_dir / "category2" / "file2.txt").touch()
        (cache_dir / ".hidden").touch()  # Should be skipped

        cleared = manager.clear_cache()

        # Should clear category dirs but not hidden files
        assert len(cleared) == 2
        assert not (cache_dir / "category1").exists()
        assert not (cache_dir / "category2").exists()
        assert (cache_dir / ".hidden").exists()  # Hidden files preserved

    def test_clear_cache_by_category(self, manager, tmp_path):
        """Test clearing specific category."""
        cache_dir = manager.cache_dir
        (cache_dir / "guidescan_indexes").mkdir()
        (cache_dir / "guidescan_indexes" / "file.txt").touch()
        (cache_dir / "other_category").mkdir()
        (cache_dir / "other_category" / "file.txt").touch()

        cleared = manager.clear_cache("guidescan_indexes")

        assert len(cleared) == 1
        assert not (cache_dir / "guidescan_indexes" / "file.txt").exists()
        assert (cache_dir / "other_category" / "file.txt").exists()

    def test_clear_cache_nonexistent_category(self, manager):
        """Test clearing non-existent category returns empty list."""
        cleared = manager.clear_cache("nonexistent")
        assert cleared == []

    def test_clear_cache_empty_directory(self, manager):
        """Test clearing empty cache directory."""
        cleared = manager.clear_cache()
        assert cleared == []


class TestListAvailableFiles:
    """Tests for list_available_files method."""

    def test_list_available_files(self, manager):
        """Test listing all available files."""
        files = manager.list_available_files()

        assert "guidescan_indexes" in files
        assert "human" in files["guidescan_indexes"]
        assert "mouse" in files["guidescan_indexes"]
        assert files["guidescan_indexes"]["human"]["organism_id"] == "NCBITaxon:9606"


class TestIntegrationWithAnnotateGuides:
    """Integration tests for use with annotate_guides module."""

    def test_fetch_returns_files_for_index_lookup(self, manager):
        """Test that fetch returns files suitable for finding .index.gs files."""
        with patch("cellxgene_schema.reference_file_manager.pooch.retrieve") as mock_retrieve:
            # Simulate extracted guidescan index files
            mock_retrieve.return_value = [
                "/cache/guidescan_indexes/human.index.forward",
                "/cache/guidescan_indexes/human.index.reverse",
                "/cache/guidescan_indexes/human.index.gs",
            ]

            files = manager.fetch("guidescan_indexes", "human")

            # Should be able to find .index.gs file
            gs_files = [f for f in files if f.endswith(".index.gs")]
            assert len(gs_files) == 1
            assert gs_files[0].endswith(".index.gs")

            # Should be able to derive index prefix
            index_prefix = gs_files[0][: -len(".index.gs")]
            assert index_prefix == "/cache/guidescan_indexes/human"

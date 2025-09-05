"""
Dataset registry for PySEE performance testing.

This module provides a registry system for managing curated datasets
with metadata, checksums, and performance targets.
"""

import yaml
import os
from pathlib import Path
from typing import Dict, Any, Optional, List
import anndata as ad
import scanpy as sc
import hashlib
import requests

# Import will be handled at runtime


class DatasetRegistry:
    """Registry for managing performance testing datasets."""

    def __init__(self, config_path: str = "configs/datasets.yml"):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.data_dir = Path("data")
        self.raw_dir = self.data_dir / "raw"
        self.processed_dir = self.data_dir / "processed"

        # Create directories
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.processed_dir.mkdir(parents=True, exist_ok=True)

    def _load_config(self) -> Dict[str, Any]:
        """Load dataset configuration from YAML file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Dataset config not found: {self.config_path}")

        with open(self.config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)

    def get_dataset_info(self, dataset_id: str) -> Dict[str, Any]:
        """Get information about a specific dataset."""
        if dataset_id not in self.config["datasets"]:
            raise KeyError(f"Dataset '{dataset_id}' not found in registry")

        return self.config["datasets"][dataset_id]

    def list_datasets(self, size_category: Optional[str] = None) -> List[str]:
        """List available datasets, optionally filtered by size category."""
        datasets = list(self.config["datasets"].keys())

        if size_category:
            datasets = [d for d in datasets if self.config["datasets"][d]["size"] == size_category]

        return datasets

    def get_datasets_by_category(self) -> Dict[str, List[str]]:
        """Get datasets grouped by size category."""
        categories = {}

        for dataset_id, info in self.config["datasets"].items():
            size = info["size"]
            if size not in categories:
                categories[size] = []
            categories[size].append(dataset_id)

        return categories

    def get_performance_targets(self, size_category: str) -> Dict[str, float]:
        """Get performance targets for a size category."""
        targets = self.config["performance_targets"]

        return {
            "memory_mb": targets["memory_usage"][size_category],
            "rendering_time": targets["rendering_time"][size_category],
            "interactive_response": targets["interactive_response"],
        }

    def validate_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """Validate file checksum."""
        if not expected_checksum or expected_checksum.startswith("sha256:"):
            # Extract hash from sha256:hash format
            expected_hash = expected_checksum.split(":", 1)[1] if ":" in expected_checksum else expected_checksum

            # Calculate file hash
            sha256_hash = hashlib.sha256()
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(chunk)

            actual_hash = sha256_hash.hexdigest()
            return actual_hash == expected_hash

        return True  # Skip validation if no checksum provided

    def download_dataset(self, dataset_id: str, force: bool = False) -> Path:
        """Download a dataset if not already present."""
        info = self.get_dataset_info(dataset_id)

        if info["source"] == "scanpy_builtin":
            # Built-in datasets don't need downloading
            return self._get_scanpy_dataset(dataset_id)

        download_url = info.get("download_url")
        if not download_url:
            raise ValueError(f"Dataset '{dataset_id}' has no download URL")

        # Determine file extension
        if download_url.endswith(".h5"):
            filename = f"{dataset_id}.h5"
        elif download_url.endswith(".h5ad"):
            filename = f"{dataset_id}.h5ad"
        else:
            filename = f"{dataset_id}.h5"  # Default to h5

        file_path = self.raw_dir / filename

        # Check if file already exists
        if file_path.exists() and not force:
            # Validate checksum if provided
            checksum = info.get("checksum")
            if checksum and not self.validate_checksum(file_path, checksum):
                print(f"⚠️ Checksum validation failed for {dataset_id}, re-downloading...")
                force = True
            else:
                print(f"✅ Dataset {dataset_id} already exists: {file_path}")
                return file_path

        # Download the dataset
        print(f"📥 Downloading {dataset_id} from {download_url}...")

        try:
            response = requests.get(download_url, stream=True)
            response.raise_for_status()

            total_size = int(response.headers.get("content-length", 0))
            downloaded = 0

            with open(file_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)

                        if total_size > 0:
                            progress = (downloaded / total_size) * 100
                            print(f"\r   Progress: {progress:.1f}%", end="", flush=True)

            print(f"\n✅ Downloaded {dataset_id}: {file_path}")

            # Validate checksum
            checksum = info.get("checksum")
            if checksum and not self.validate_checksum(file_path, checksum):
                print(f"⚠️ Warning: Checksum validation failed for {dataset_id}")

            return file_path

        except Exception as e:
            print(f"❌ Failed to download {dataset_id}: {e}")
            if file_path.exists():
                file_path.unlink()  # Remove partial download
            raise

    def _get_scanpy_dataset(self, dataset_id: str) -> Path:
        """Get scanpy built-in dataset and save to cache."""
        file_path = self.processed_dir / f"{dataset_id}.h5ad"

        if file_path.exists():
            return file_path

        print(f"📥 Loading scanpy dataset: {dataset_id}")

        if dataset_id == "pbmc3k":
            adata = sc.datasets.pbmc3k()
        elif dataset_id == "pbmc68k":
            adata = sc.datasets.pbmc68k_reduced()
        else:
            raise ValueError(f"Unknown scanpy dataset: {dataset_id}")

        # Save to cache
        adata.write(file_path)
        print(f"✅ Cached {dataset_id}: {file_path}")

        return file_path

    def load_dataset(self, dataset_id: str, backed: bool = False, check_requirements: bool = True) -> ad.AnnData:
        """Load a dataset from cache or download if needed."""
        info = self.get_dataset_info(dataset_id)

        # Check system requirements if requested
        if check_requirements:
            self._check_system_requirements(dataset_id, info)

        # Check if it's a synthetic dataset
        if info["source"] == "synthetic":
            return self._get_synthetic_dataset(dataset_id)

        # Check if it's a scanpy built-in dataset
        if info["source"] == "scanpy_builtin":
            file_path = self._get_scanpy_dataset(dataset_id)
        else:
            # Download if needed
            file_path = self.download_dataset(dataset_id)

        # Load the dataset
        print(f"📊 Loading dataset: {dataset_id}")

        if backed and info["size"] in ["large", "very_large"]:
            # Use backed mode for large datasets
            adata = ad.read_h5ad(file_path, backed="r")
            print(f"✅ Loaded {dataset_id} in backed mode")
        else:
            adata = ad.read_h5ad(file_path)
            print(f"✅ Loaded {dataset_id} in memory")

        return adata

    def _check_system_requirements(self, dataset_id: str, info: Dict[str, Any]):
        """Check system requirements for a dataset."""
        try:
            # Import here to avoid circular imports
            import sys
            from pathlib import Path

            project_root = Path(__file__).parent.parent.parent
            sys.path.insert(0, str(project_root))

            from pysee.utils.system_requirements import check_system_requirements

            # Check requirements
            compatible = check_system_requirements(dataset_id, info["size"], info["memory_mb"])

            if not compatible:
                print(f"⚠️ Proceeding with {dataset_id} despite memory warnings...")

        except ImportError:
            # System requirements checker not available
            pass
        except Exception as e:
            print(f"⚠️ Could not check system requirements: {e}")

    def _get_synthetic_dataset(self, dataset_id: str) -> ad.AnnData:
        """Get synthetic dataset."""
        # Import here to avoid circular imports
        from .dataset_fixtures import DatasetFixtures

        if dataset_id == "synthetic_small":
            return DatasetFixtures.generate_synthetic_small()
        elif dataset_id == "synthetic_medium":
            return DatasetFixtures.generate_synthetic_medium()
        elif dataset_id == "synthetic_large":
            return DatasetFixtures.generate_synthetic_large()
        elif dataset_id == "synthetic_very_large":
            return DatasetFixtures.generate_synthetic_very_large()
        else:
            raise ValueError(f"Unknown synthetic dataset: {dataset_id}")

    def get_all_datasets(self, size_categories: Optional[List[str]] = None) -> Dict[str, ad.AnnData]:
        """Get all datasets, optionally filtered by size categories."""
        datasets = {}

        if size_categories:
            dataset_ids = []
            for category in size_categories:
                dataset_ids.extend(self.list_datasets(category))
        else:
            dataset_ids = self.list_datasets()

        for dataset_id in dataset_ids:
            try:
                datasets[dataset_id] = self.load_dataset(dataset_id)
            except Exception as e:
                print(f"⚠️ Failed to load {dataset_id}: {e}")
                continue

        return datasets

    def print_summary(self):
        """Print a summary of available datasets."""
        print("📊 PySEE Dataset Registry Summary")
        print("=" * 50)

        categories = self.get_datasets_by_category()

        for category, dataset_ids in categories.items():
            print(f"\n{category.upper()} Datasets:")
            for dataset_id in dataset_ids:
                info = self.get_dataset_info(dataset_id)
                print(f"  {dataset_id:20} | {info['cells']:8,} cells | {info['genes']:6,} genes | {info['memory_mb']:6.0f} MB")

        print(f"\nTotal datasets: {len(self.config['datasets'])}")
        print(f"Data directory: {self.data_dir.absolute()}")


def main():
    """Test the dataset registry."""
    registry = DatasetRegistry()
    registry.print_summary()

    # Test loading a small dataset
    print("\n🧪 Testing dataset loading...")
    try:
        adata = registry.load_dataset("pbmc3k")
        print(f"✅ Successfully loaded pbmc3k: {adata.n_obs} cells, {adata.n_vars} genes")
    except Exception as e:
        print(f"❌ Failed to load pbmc3k: {e}")


if __name__ == "__main__":
    main()

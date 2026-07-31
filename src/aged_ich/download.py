from __future__ import annotations

from pathlib import Path
from typing import Any

import requests

from .config import (
    AGE_EFFECT_DATA_DIRECTORY,
    AGE_EFFECT_GITHUB_REPOSITORIES,
)


def _github_contents_url(repository: str) -> str:
    return (
        f"https://api.github.com/repos/{repository}/contents/"
        f"{AGE_EFFECT_DATA_DIRECTORY}"
    )


def list_cns_files(
    session: requests.Session,
    repository: str,
    timeout_seconds: int = 60,
) -> list[dict[str, Any]]:
    """List the CNS age-regression tables available in the source repository."""

    response = session.get(
        _github_contents_url(repository),
        timeout=timeout_seconds,
    )
    response.raise_for_status()

    items = response.json()
    files = [
        item
        for item in items
        if item.get("type") == "file"
        and item.get("name", "").startswith("brain_")
        and item.get("name", "").endswith("_pvalue_results.csv")
    ]

    return sorted(files, key=lambda item: item["name"])


def download_cns_age_tables(
    output_dir: Path,
    repositories: tuple[str, ...] = AGE_EFFECT_GITHUB_REPOSITORIES,
    timeout_seconds: int = 60,
) -> list[Path]:
    """
    Download all brain/CNS age-regression CSV files from Age Effect Explorer.

    The function tries the current and historical repository names. It expects
    13 CNS files and returns their local paths.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update(
        {"User-Agent": "aged-ich-target-prioritization/1.0"}
    )

    errors: list[str] = []

    for repository in repositories:
        try:
            file_metadata = list_cns_files(
                session=session,
                repository=repository,
                timeout_seconds=timeout_seconds,
            )

            if len(file_metadata) != 13:
                raise RuntimeError(
                    f"Expected 13 CNS tables in {repository}; found {len(file_metadata)}."
                )

            downloaded: list[Path] = []
            for item in file_metadata:
                destination = output_dir / item["name"]

                if not destination.exists():
                    response = session.get(
                        item["download_url"],
                        timeout=timeout_seconds,
                    )
                    response.raise_for_status()
                    destination.write_bytes(response.content)

                downloaded.append(destination)

            return downloaded

        except (requests.RequestException, RuntimeError, KeyError, TypeError) as error:
            errors.append(f"{repository}: {error}")

    error_text = "\n".join(errors)
    raise RuntimeError(
        "Could not download the Age Effect Explorer CNS tables.\n"
        f"{error_text}"
    )

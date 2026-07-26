"""
Internal HTTP helper shared by every module in the package.

Centralizes the ISDA base URL, request/timeout handling, and logging so
individual modules don't repeat the same try/except boilerplate.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

import requests

logger = logging.getLogger("pyisda")
if not logger.handlers:
    # Library-friendly default: don't spam stdout unless the caller
    # configures logging themselves, but make it easy to opt in.
    logger.addHandler(logging.NullHandler())

#: Base URL for the IBDC ISDA REST API. IBDC migrated from
#: ibdc.dbtindia.gov.in to ibdc.dbt.gov.in; this points at the current
#: domain. Can be overridden at runtime, e.g. for a staging environment:
#:     import pyisda
#:     pyisda.config.ISDA_BASE_URL = "https://staging.example.org/isda/api"
ISDA_BASE_URL = "https://ibdc.dbt.gov.in/isda/api"

DEFAULT_TIMEOUT = 10


class ISDARequestError(Exception):
    """Raised when a request to the ISDA API fails outright (network/HTTP)."""


def get_json(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    GET a URL and return the parsed JSON body.

    Raises:
        ISDARequestError: on timeout, connection error, HTTP error status,
            or a response body that isn't valid JSON. Callers in this
            package catch this and convert it into a logged warning plus
            a ``None``/empty return, matching the rest of the API.
    """
    try:
        response = requests.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.Timeout as exc:
        raise ISDARequestError(f"Request timed out: {url}") from exc
    except requests.exceptions.ConnectionError as exc:
        raise ISDARequestError(f"Connection error for URL: {url}") from exc
    except requests.exceptions.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else "?"
        reason = exc.response.reason if exc.response is not None else ""
        raise ISDARequestError(f"HTTP error {status} {reason} for URL: {url}") from exc
    except ValueError as exc:  # json.JSONDecodeError subclasses ValueError
        raise ISDARequestError(f"Failed to decode JSON response from: {url}") from exc


def download_file(url: str, dest_path: str, timeout: int = 30) -> str:
    """
    Download a URL's raw body (e.g. a structure file) to `dest_path`.

    Args:
        url: URL to fetch.
        dest_path: Local file path to write the response body to.
        timeout: Request timeout in seconds.

    Returns:
        `dest_path`, for convenient chaining.

    Raises:
        ISDARequestError: on timeout, connection error, HTTP error status,
            or an empty response body.
    """
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        if not response.content:
            raise ISDARequestError(f"Empty response body from: {url}")
        with open(dest_path, "wb") as fh:
            fh.write(response.content)
        return dest_path
    except requests.exceptions.Timeout as exc:
        raise ISDARequestError(f"Request timed out: {url}") from exc
    except requests.exceptions.ConnectionError as exc:
        raise ISDARequestError(f"Connection error for URL: {url}") from exc
    except requests.exceptions.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else "?"
        reason = exc.response.reason if exc.response is not None else ""
        raise ISDARequestError(f"HTTP error {status} {reason} for URL: {url}") from exc

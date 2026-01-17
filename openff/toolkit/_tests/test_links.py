import pathlib
import re
from urllib.request import Request, urlopen

import pytest


def find_readme_links() -> list[str]:
    """Yield all the links in the main README.md file.

    Returns
    -------
    readme_examples : list[str]
        The list of links included in the README.md file.
    """
    if "site-packages" in __file__:
        # This test file is being collected from the installed package, which
        # does not provide the README file.
        # Note that there will likely be a mis-bundled file
        # $CONDA_PREFIX/lib/python3.x/site-packages/README.md, but this is not
        # the toolkit's README file!
        return list()

    else:
        readme_file_path = pathlib.Path(__file__).parents[3] / "README.md"
        with open(readme_file_path.as_posix()) as f:
            readme_content = f.read()

    with open(readme_file_path.as_posix()) as f:
        readme_content = f.read()

    return re.findall("http[s]?://(?:[0-9a-zA-Z]|[-/.%:_])+", readme_content)


@pytest.mark.parametrize("readme_link", find_readme_links())
def test_readme_links(readme_link):
    """Test URLs in README are reachable"""
    # Some websites do not accept requests that don't specify the
    # client and the type of accepted documents so we add fake info
    # to avoid the response being an error.
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/xhtml+xml,text/html,application/xml;q=0.9,*/*;q=0.8",
    }
    request = Request(readme_link, headers=headers)

    # Some DOI-based links are now behind DDoS protections, so skip them
    if "doi.org" in readme_link:
        pytest.skip("DOI links are behind DDoS protection and do not resolve")
    if "codecov.io" in readme_link:
        pytest.skip("Codecov website DOI also may be behind DDoS protection")
    if "zenodo" in readme_link:
        # October 2023 Zenodo upgrade seems to have broken the bibtex export
        # link although it is still accessible from the web UI (go to the
        # landing page of a record, Export section (select BibTeX in
        # drop-down) and click "Export", will work in browser)

        # Jan 2026 Zenodo seems to have blocked programmatic (non-browser)
        # access to regular release pages as well, so skip all zenodo checks.
        # https://support.zenodo.org/help/en-gb/13-policies/29-why-is-my-ip-address-blocked
        pytest.skip("Programmatic access tests no longer work for Zenodo URLs")
    if readme_link.endswith("MIT"):
        # January 22 2024: Link causing some failures, accessible in browser
        pytest.skip("Flaky")

    # Try to connect 5 times, keeping track of exceptions so useful feedback can be provided.
    success = False
    exception = None
    for retry in range(5):
        try:
            urlopen(request)
            success = True
            break
        except Exception as e:
            exception = e
    if not (success):
        raise exception

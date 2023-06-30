import os
import re
from urllib.request import Request, urlopen

import pytest

ROOT_DIR_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", ".."
)


def find_readme_links():
    """Yield all the links in the main README.md file.

    Returns
    -------
    readme_examples : List[str]
        The list of links included in the README.md file.
    """
    readme_file_path = os.path.join(ROOT_DIR_PATH, "README.md")
    with open(readme_file_path, "r") as f:
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

    # Try to connect 5 times, keeping track of exceptions so useful feedback can be provided.
    success = False
    exception = None
    for retry in range(5):  # noqa: B007
        try:
            urlopen(request)
            success = True
            break
        except Exception as e:
            exception = e
    if not (success):
        raise exception

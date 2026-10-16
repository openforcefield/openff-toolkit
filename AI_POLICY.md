# AI Policy

## OpenFF Ecosystem Policies

The Open Force Field Initiative reserves the right to reject contributions if we judge that they do not follow the spirit of these principles.
The following policies hold for all repositories in the OpenFF ecosystem:

1. Repository Maintainers always reserve the right to make exceptions to any of the following rules on a case-by-case basis.
2. LLMs use must be always disclosed. At minimum in the PR description, but preferably reflected in a “Co-Authored-By” line in the commit.
3. LLMs may always be used for translation (disclosing LLM use is important here, so that maintainers can be more considerate when encountering translation inaccuracies).
4. All content must have a human identified as being accountable for the contribution. The human takes responsibility for the contribution and stands by its correctness. No agent-only commits.
5. Code authors must retain full authorship for any contributed code and contributions must be covered by the [MIT license](https://github.com/openforcefield/openff-toolkit/blob/main/LICENSE), and [Developer Certificate of Origin](https://developercertificate.org/).

## Policies specific to this repository

The following policies are specific only to this repository, and decided by the repository's maintainer(s).

In this repository, it is either acceptable (✅) or not acceptable (❌) to use LLMs for:

- ✅️ Generating content that only the author sees, but affects their work (searching repositories, asking questions, finding bugs, self-review, linting).
- ✅️ Generating code that is shared, but not committed to the code base (e.g. Minimally Reproducible Examples in issue text).
- Generating content which the human author has:

|                                                              | Code | Docs | PRs & Issues | Replies to humans |
|--------------------------------------------------------------|:----:|:----:|:------------:|:-----------------:|
| Read thoroughly and understands every line before submission |  ✅️  |  ✅️  |      ✅️      |         ❌         |
| Checked for obvious issues before submission                 |  ❌   |  ❌   |      ❌       |         ❌         |
| Not reviewed before submission                               |  ❌   |  ❌   |      ❌       |         ❌         |

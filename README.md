# RSeQC

[![PyPI version](https://img.shields.io/pypi/v/RSeQC.svg)](https://pypi.org/project/RSeQC/)
[![Python Versions](https://img.shields.io/pypi/pyversions/RSeQC.svg)](https://pypi.org/project/RSeQC/)
[![License](https://img.shields.io/github/license/liguowang/RSeQC)](LICENSE)
[![CI](https://github.com/liguowang/RSeQC/actions/workflows/ci.yml/badge.svg)](https://github.com/liguowang/RSeQC/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://www.lwang.org/RSeQC/)
[![PyPI Downloads](https://static.pepy.tech/badge/rseqc)](https://pepy.tech/projects/rseqc)

## Requirements

- Python 3.9 or later
- `pip`

Check your Python and pip installations:

```bash
python3 --version
python3 -m pip --version
```

Example:

```text
Python 3.11.13
pip 25.3 from /path/to/lib/python3.11/site-packages/pip (python 3.11)
```

If `pip` is not installed, follow the official installation instructions:

https://pip.pypa.io/en/stable/installation/

---

## (Optional) Create a Virtual Environment

Using a virtual environment is recommended to keep RSeQC and its dependencies isolated from other Python packages.

Create a new environment:

```bash
python3 -m venv rseqc_env
```

Activate it:

**Linux/macOS**

```bash
source rseqc_env/bin/activate
```

**Windows (PowerShell)**

```powershell
rseqc_env\Scripts\Activate.ps1
```

---

## Install RSeQC

### Install the latest stable release from [PyPI](https://pypi.org/project/RSeQC/)

```bash
python3 -m pip install rseqc
```

### Install the latest development version from GitHub

```bash
python3 -m pip install git+https://github.com/liguowang/RSeQC.git
```

You can also install a specific branch, tag, or commit.

Install a specific branch:

```bash
python3 -m pip install git+https://github.com/liguowang/RSeQC.git@main
```

Install a specific release tag:

```bash
python3 -m pip install git+https://github.com/liguowang/RSeQC.git@v5.0.6
```

---

## Upgrade RSeQC

Upgrade an existing installation from PyPI:

```bash
python3 -m pip install --upgrade rseqc
```

To upgrade a GitHub installation:

```bash
python3 -m pip install --upgrade git+https://github.com/liguowang/RSeQC.git
```

---

## Uninstall RSeQC

Remove RSeQC:

```bash
python3 -m pip uninstall rseqc
```

---

## Verify the Installation

Check that RSeQC is installed correctly:

```bash
bam_stat.py --help
```

or

```bash
python3 -c "from importlib.metadata import version; print('RSeQC', version('RSeQC'))"
```

---

## License

RSeQC is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

See the `LICENSE` file included with this project for the complete license text.

---

## Documentation

Documentation and tutorials are available at:

1. http://rseqc.sourceforge.net/
2. https://www.lwang.org/RSeQC/

## References

1. Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012;28(16):2184-2185. doi:10.1093/bioinformatics/bts356
2. Wang L, Nie J, Sicotte H, et al. Measure transcript integrity using RNA-seq data. BMC Bioinformatics. 2016;17:58. doi:10.1186/s12859-016-0922-z


# RSeQC

## Requirements

- Python 3.9 or later
- `pip`

Check your Python and pip installations:

```bash
python3 --version
pip --version
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

Creating a virtual environment is recommended to avoid conflicts with other Python packages.

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

### Install the latest stable release from PyPI

```bash
pip install rseqc
```

### Install the latest development version from GitHub

```bash
pip install git+https://github.com/liguowang/RSeQC.git
```

You can also install a specific branch, tag, or commit.

Install a specific branch:

```bash
pip install git+https://github.com/liguowang/RSeQC.git@main
```

Install a specific release tag:

```bash
pip install git+https://github.com/liguowang/RSeQC.git@v5.05
```

---

## Upgrade RSeQC

Upgrade an existing installation from PyPI:

```bash
pip install --upgrade rseqc
```

To upgrade a GitHub installation:

```bash
pip install --upgrade git+https://github.com/liguowang/RSeQC.git
```

---

## Uninstall RSeQC

Remove RSeQC:

```bash
pip uninstall rseqc
```

---

## Verify the Installation

Check that RSeQC is installed correctly:

```bash
bam_stat.py --help
```

or

```bash
python -c "import qcmodule; print('RSeQC installed successfully')"
```

---

## License

RSeQC is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

See the `LICENSE` file included with this project for the complete license text.

---

# Documentation

Documentation and tutorials are available at:

http://rseqc.sourceforge.net/

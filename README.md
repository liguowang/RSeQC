# Installing RSeQC

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

Install the latest release from PyPI:

```bash
pip install rseqc
```

---

## Upgrade RSeQC

Upgrade to the latest version:

```bash
pip install --upgrade rseqc
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

RSeQC is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

---

# Documentation

Documentation and tutorials are available at:

http://rseqc.sourceforge.net/

# Install RSeQC using [pip](https://pip.pypa.io/en/stable/)

### Check if [pip](https://pip.pypa.io/en/stable/) is available

`pip` is a package manager that allows users to install and manage Python packages. It comes with Python 3.4 and later versions.
 
 `pip --version`
 
 `pip 24.2 from /XXX/YYY/lib/python3.11/site-packages/pip (python 3.11)`

 If you see any errors, following this [instructions](https://pypi.org/project/pip/) to install `pip`.

### (Optional) Create Virtual Environments (Note: `venv` is available in Python 3.3 and later. You can also use [virtualenv](https://packaging.python.org/en/latest/key_projects/#virtualenv))

 `$ python3 -m venv my_env` (will create a directory called my_env)

 `$ source my_env/bin/activate`

### Install

 `$ pip install rseqc`

### Upgrade

 `$ pip install rseqc --upgrade`

### Uninstall

 `pip -y uninstall rseqc`

# Documentation

http://rseqc.sourceforge.net/

# PGS Grep

PGS Grep is a utility for searching SNPs within Polygenic Score (PGS) files from the PGS Catalog (https://www.pgscatalog.org/). It helps users locate variant information across PGS datasets and optionally related results using linkage equilibrium (LD) data from the ([1000 Genome project](https://www.internationalgenome.org/)). This LD information is fetched from the public ([LDlink API](https://ldlink.nih.gov/apiaccess)).

## Getting Started

These instructions will help you set up the project and run the search tool on your local machine.

### Prerequisites

* Python 3.8 or newer
* An API key from LDlink: https://ldlink.nih.gov/apiaccess
* All dependencies listed in `requirements.txt`

### Installing (This needs to be revised)

1. Clone the repository.
2. Create and activate a Python virtual environment.
3. Install dependencies from `requirements.txt`.

```
python -m venv venv
source venv/bin/activate   # or venv\Scripts\activate on Windows
pip install -r requirements.txt
```

## Usage (add more details, edit so that command can run that launches startup)
First, find and download the PGS file you want to scan from the PGS Catalog https://www.pgscatalog.org/


Run this command to create the UI window:

```
if ((Split-Path $PWD -Leaf) -ne "src") { cd src }
python -m streamlit run StreamlitGUI.py
```

The script scans a PGS file for matching SNP entries and reports relevant coordinates and metadata. When LDlink token access is configured, it can also retrieve LD information for selected variants.



## Running the tests

To run tests, run 
```
pip install pytest pytest-mock requests
python -m pytest StreamlitTestSuite.py -v
```

## Built With

* Python 3

## Contributing

Contributions are welcome. Submit bug reports and feature requests through the repository issue tracker.

## Authors

* **Project Maintainer** - *Current work*

## License (This needs to be revised)

PGS Grep is released under the GNU General Public License v3.0 (GPL-3.0). You are free to reuse and modify this code under the terms of the GPL-3.0 license, but any derivative work must also be released under the same license. See details at https://www.gnu.org/licenses/gpl-3.0.html

## Acknowledgments

* PGS Catalog for sourcing polygenic score data
* LDlink for LD data access
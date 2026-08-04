# PGS Grep

PGS Grep is a utility for searching SNPs within Polygenic Score (PGS) files from the [PGS Catalog](https://www.pgscatalog.org/). It helps users locate variant information across PGS datasets and, optionally, related results using linkage disequilibrium (LD) data from the [1000 Genomes Project](https://www.internationalgenome.org/). This LD information is fetched from the public [LDlink API](https://ldlink.nih.gov/apiaccess).

## Getting Started

These instructions will help you set up the project and run the search tool on your local machine.

### Prerequisites

* **Python 3.10 or newer**
* All dependencies listed in `requirements.txt`
* *Optional:* a free personal token from LDlink (https://ldlink.nih.gov/apiaccess), needed only if you want to search LD proxy variants. The core scan works without one.

### Installing

```
git clone https://github.com/CharlieBuhanan/PGS-Grep-Public.git
cd PGS-Grep-Public

python -m venv .venv
.venv\Scripts\activate          # Windows
# source .venv/bin/activate     # macOS / Linux

pip install -r requirements.txt
```

## Usage

First, find and download the score you want to scan from the PGS Catalog (https://www.pgscatalog.org/). Download the **harmonized** scoring file — the one whose filename ends in `_hmPOS_GRCh38.txt.gz` or `_hmPOS_GRCh37.txt.gz`, depending on which genome build you want to search. Leave it gzipped; the app reads `.gz` directly and accepts uploads up to 200 MB. Optionally download the matching `.md5` checksum file as well, and the app will verify the download's integrity for you.

Then, **from the repository root**, launch the app:

```
python -m streamlit run StreamlitGUI.py
```

Running from the repository root matters: Streamlit loads `.streamlit/config.toml` and serves the bundled Roboto fonts from `static/` relative to the current working directory, so launching from elsewhere will lose the app's theming.

This opens a guided wizard in your browser that walks through uploading the scoring file, entering the rsID or chromosome/position you want to find, optionally supplying an LDlink token to include LD proxies, and running the scan. Results can be reviewed in the browser and exported as CSV.

## Running the tests

```
pip install pytest
python -m pytest StreamlitTestSuite.py -v
```

## Built With

* [Python 3](https://www.python.org/) (3.10+)
* [Streamlit](https://streamlit.io/) — the wizard user interface
* [pandas](https://pandas.pydata.org/) — scoring file parsing and results tables
* [Requests](https://requests.readthedocs.io/) — LDlink API access

## Contributing

Contributions are welcome. Submit bug reports and feature requests through the [issue tracker](https://github.com/CharlieBuhanan/PGS-Grep-Public/issues).

## Authors

* **Charlie W. Buhanan** — *Author*
* **Deborah H. Glueck** — *Supervision*

## License

PGS Grep is released under the GNU General Public License v3.0 (GPL-3.0). You are free to reuse and modify this code under the terms of the GPL-3.0 license, but any derivative work must also be released under the same license. See the [LICENSE](LICENSE) or https://www.gnu.org/licenses/gpl-3.0.html for details.

## Acknowledgments

* [PGS Catalog](https://www.pgscatalog.org/) for sourcing polygenic score data
* [LDlink](https://ldlink.nih.gov/), from the National Cancer Institute, for LD data access
* [1000 Genomes Project](https://www.internationalgenome.org/) for the underlying reference population data
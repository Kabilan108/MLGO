# --------------------------------- VARIABLES ---------------------------------
RSCRIPT = Rscript
PYTHON = python


# -------------------------------- ENV SETUP ----------------------------------
git_init:
	git init
	git config --global user.name "Muaddib"
	git config --global user.email "tonykabilanokeke@gmail.com"
	git config --global core.editor "nvim"

install: 
	@echo "Installing..."
	poetry install
	poetry run pre-commit install

setup: git_init install


# --------------------------------- TESTING -----------------------------------
test:
	pytest


# --------------------------------- DOCS --------------------------------------
docs_view:
	@echo View API documentation... 
	PYTHONPATH=src pdoc src --http localhost:8080

docs_save:
	@echo Save documentation to docs... 
	PYTHONPATH=src pdoc src -o docs


# ----------------------------- DATA PIPELINE ---------------------------------
data/raw/DEG.rds: scripts/loaddatasets.R
	@echo "Loading data..."
	$(RSCRIPT) scripts/loaddatasets.R

data/processed/clean-data.feather: data/raw/DEG.rds scripts/clean_data.R
	@echo "Cleaning data..."
	$(RSCRIPT) scripts/clean_data.R

data/final/xy.pkl: data/processed/clean-data.feather scripts/prepdata.py
	@echo "Preparing ML data..."
	$(PYTHON) scripts/prepdata.py

get_data: data/final/xy.pkl

pipeline: data/processed/xy.pkl models/svc.pkl notebooks/results.ipynb


# ------------------------------- CLEAN UP ------------------------------------
.PHONY: clean

clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete
	rm -rf .pytest_cache

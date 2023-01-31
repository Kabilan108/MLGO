# Environment Setup
git_init:
	git init
	git config --global user.name "Muaddib"
	git config --global user.email "tonykabilanokeke@gmail.com"
	git config --global core.editor "nvim"

install: 
	@echo "Installing..."
	poetry install
	poetry run pre-commit install

download_data:
	@echo "Downloading data..."
	wget https://gist.githubusercontent.com/khuyentran1401/a1abde0a7d27d31c7dd08f34a2c29d8f/raw/da2b0f2c9743e102b9dfa6cd75e94708d01640c9/Iris.csv -O data/raw/iris.csv

setup: git_init install download_data


# Testing
test:
	pytest


# Documentation
docs_view:
	@echo View API documentation... 
	PYTHONPATH=src pdoc src --http localhost:8080

docs_save:
	@echo Save documentation to docs... 
	PYTHONPATH=src pdoc src -o docs


# Data processing
data/processed/xl.pkl: data/raw src/process.py
	@echo "Processing data..."
	python src/process.py

models/svc.pkl: data/processed/xl.pkl src/train_model.py
	@echo "Training model..."
	python src/train_model.py

pipeline: data/processed/xl.pkl models/svc.pkl


# data/processed/xy.pkl: data/raw src/process.py
# 	@echo "Processing data..."
# 	python src/process.py
#
# models/svc.pkl: data/processed/xy.pkl src/train_model.py
# 	@echo "Training model..."
# 	python src/train_model.py
#
# notebooks/results.ipynb: models/svc.pkl src/run_notebook.py
# 	@echo "Running notebook..."
# 	python src/run_notebook.py

# pipeline: data/processed/xy.pkl models/svc.pkl notebooks/results.ipynb


# Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete
	rm -rf .pytest_cache


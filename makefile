.PHONY: init env docs live-docs clean

PIP := $(shell command -v uv > /dev/null && echo "uv pip" || echo "pip")

init:
	mkdir -p data models
	cp .env.template .env
	$(CONDA_EXE) env create -f env.yml
	$(CONDA_EXE) run -n {{cookiecutter.short_name}} $(PIP) install -r requirements.txt
	$(CONDA_EXE) run -n {{cookiecutter.short_name}} conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y
	$(CONDA_EXE) run -n {{cookiecutter.short_name}} $(PIP) install -e .
	$(CONDA_EXE) run -n {{cookiecutter.short_name}} pre-commit install

env:
	$(PIP) freeze > requirements.txt

docs:
	quarto render docs/index.qmd --to html

live-docs:
	quarto preview docs/blog-post.qmd --to html --port 7777 --no-browser

clean:
	rm -rf models/*
	rm -rf output/*
	rm -rf wandb/*
	rm -rf docs/.quarto docs/*.html
	rm -rf src/*.egg-info

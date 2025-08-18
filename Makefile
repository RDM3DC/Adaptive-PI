PYTHON ?= python3
VENV ?= .venv
PIP := $(VENV)/bin/pip
PY := $(VENV)/bin/python

.PHONY: help venv install test lint clean phase8 phase9 phase10 all

help:
	@echo 'Targets:'
	@echo '  venv      - create virtual environment'
	@echo '  install   - install requirements'
	@echo '  test      - run pytest'
	@echo '  phase8    - run phase 8 script'
	@echo '  phase9    - run phase 9 script'
	@echo '  phase10   - run phase 10 script'
	@echo '  clean     - remove build artifacts'

venv:
	@test -d $(VENV) || $(PYTHON) -m venv $(VENV)
	@echo 'Virtualenv ready.'

install: venv
	$(PIP) install --upgrade pip
	@if [ -f requirements.txt ]; then $(PIP) install -r requirements.txt; fi
	@echo 'Dependencies installed.'

test: install
	$(PY) -m pytest -q

phase8: install
	bash scripts/run_phase8_m1.sh

phase9: install
	bash scripts/run_phase9_m1.sh

phase10: install
	bash scripts/run_phase10_m1.sh

all: test

lint:
	$(PIP) install ruff
	$(VENV)/bin/ruff check src

clean:
	rm -rf $(VENV) build dist *.egg-info __pycache__ .pytest_cache
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	@echo 'Cleaned.'
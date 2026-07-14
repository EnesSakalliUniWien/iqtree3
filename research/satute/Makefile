SHELL := /bin/bash
PROJECT_ROOT := $(abspath .)
REPOSITORY_ROOT := $(abspath ../..)
PYTHON ?= $(if $(PYTHON_BIN),$(PYTHON_BIN),$(shell for candidate in python3 "$$HOME/miniconda3/bin/python3" /usr/bin/python3; do "$$candidate" -c 'import numpy, scipy, matplotlib' >/dev/null 2>&1 && { command -v "$$candidate" 2>/dev/null || echo "$$candidate"; break; }; done))
IQTREE ?= $(REPOSITORY_ROOT)/build/iqtree3
EXP ?= 003_invariant_mixture
PROFILE ?= smoke

export PYTHONPATH := $(PROJECT_ROOT)/src:$(PYTHONPATH)

.PHONY: help status audit smoke check experiment figures manuscript release clean

help:
	@echo "SatuTe research targets"
	@echo "  make status"
	@echo "  make audit"
	@echo "  make smoke"
	@echo "  make check"
	@echo "  make experiment EXP=003_invariant_mixture PROFILE=full"
	@echo "  make figures EXP=003_invariant_mixture"
	@echo "  make manuscript"
	@echo "  make release EXP=003_invariant_mixture VERSION=v1"

status:
	$(PYTHON) tools/maintenance/status.py

audit:
	$(PYTHON) tools/maintenance/audit.py

smoke:
	$(PYTHON) tools/maintenance/run_experiment.py --experiment 001_paper_reconstruction --profile smoke
	$(PYTHON) tools/maintenance/run_experiment.py --experiment 002_relative_weighting --profile smoke
	$(PYTHON) tools/maintenance/run_experiment.py --experiment 003_invariant_mixture --profile smoke
	$(PYTHON) tools/maintenance/run_experiment.py --experiment 004_saturation_diagnostics --profile smoke

check:
	$(PYTHON) tests/unit/test_head_to_head_calculations.py
	$(PYTHON) tests/unit/test_multiple_testing.py
	$(PYTHON) tests/exact/test_weighted_likelihood.py
	$(PYTHON) tests/exact/test_saturation_null.py --replicates 2000 --output artifacts/local/check/exact-null.tsv
	PYTHON_BIN=$(PYTHON) tests/native/run_satute_phase1_checks.sh $(IQTREE) artifacts/local/check/native
	$(MAKE) audit PYTHON=$(PYTHON)

experiment:
	$(PYTHON) tools/maintenance/run_experiment.py --experiment $(EXP) --profile $(PROFILE)

figures:
	$(PYTHON) tools/maintenance/run_experiment.py --experiment $(EXP) --profile full --figures

manuscript:
	mkdir -p manuscript/build
	cd manuscript && latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
	cp manuscript/build/main.pdf manuscript/manuscript.pdf

release:
	$(PYTHON) tools/maintenance/create_release_manifest.py --experiment $(EXP) --version $(VERSION)

clean:
	rm -rf manuscript/build/* artifacts/local/smoke .pytest_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

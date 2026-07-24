SHELL := /bin/bash
PROJECT_ROOT := $(abspath .)
REPOSITORY_ROOT := $(abspath ../..)
PYTHON ?= $(if $(PYTHON_BIN),$(PYTHON_BIN),$(shell for candidate in python3 "$$HOME/miniconda3/bin/python3" /usr/bin/python3; do "$$candidate" -c 'import numpy, scipy, matplotlib' >/dev/null 2>&1 && { command -v "$$candidate" 2>/dev/null || echo "$$candidate"; break; }; done))
IQTREE ?= $(REPOSITORY_ROOT)/build/iqtree3
EXP ?= 003_invariant_mixture
PROFILE ?= smoke
STUDY ?=

export PYTHONPATH := $(PROJECT_ROOT)/src:$(PYTHONPATH)

.PHONY: help status map audit plot-check smoke check power-check rate-category-check data experiment figures manuscript release clean

help:
	@echo "SatuTe research targets"
	@echo "  make status"
	@echo "  make map"
	@echo "  make audit"
	@echo "  make plot-check"
	@echo "  make smoke"
	@echo "  make check"
	@echo "  make power-check"
	@echo "  make rate-category-check"
	@echo "  make data PROFILE=smoke"
	@echo "  make publication-figures"
	@echo "  make experiment EXP=002_relative_weighting STUDY=gtr_power PROFILE=full"
	@echo "  make figures EXP=003_invariant_mixture"
	@echo "  make manuscript"
	@echo "  make release EXP=003_invariant_mixture VERSION=v1"

status:
	$(PYTHON) tools/maintenance/status.py

map:
	$(PYTHON) generate_experiment_data.py --map

audit:
	$(PYTHON) tools/maintenance/audit.py

plot-check:
	$(PYTHON) tools/maintenance/audit_power_plots.py

smoke:
	$(PYTHON) generate_experiment_data.py --profile smoke

check:
	$(PYTHON) tests/unit/test_experiment_catalog.py
	$(PYTHON) tests/unit/test_head_to_head_calculations.py
	$(PYTHON) tests/unit/test_benchmark_contracts.py
	$(PYTHON) tests/unit/test_multiple_testing.py
	$(PYTHON) tests/unit/test_fdr_calibration.py
	$(PYTHON) tests/unit/test_fdr_postprocess.py
	$(PYTHON) tests/exact/test_weighted_likelihood.py
	$(PYTHON) tests/exact/test_saturation_null.py --replicates 2000 --output artifacts/local/check/exact-null.tsv
	PYTHON_BIN=$(PYTHON) tests/native/run_satute_phase1_checks.sh $(IQTREE) artifacts/local/check/native
	$(MAKE) audit PYTHON=$(PYTHON)

power-check:
	$(PYTHON) tests/unit/test_head_to_head_calculations.py
	$(PYTHON) tests/unit/test_benchmark_contracts.py
	$(PYTHON) tests/unit/test_multiple_testing.py
	$(PYTHON) tests/unit/test_fdr_calibration.py
	$(PYTHON) tests/unit/test_fdr_postprocess.py
	$(PYTHON) tests/exact/test_weighted_likelihood.py
	$(MAKE) rate-category-check PYTHON=$(PYTHON)

rate-category-check:
	$(PYTHON) tests/native/test_satute_rate_categories.py $(IQTREE) artifacts/local/check/rate-categories

data:
	$(PYTHON) generate_experiment_data.py --profile $(PROFILE)

publication-figures:
	$(PYTHON) generate_experiment_data.py --figures

experiment:
	$(PYTHON) tools/maintenance/run_experiment.py --experiment $(EXP) --profile $(PROFILE) $(if $(STUDY),--study $(STUDY),)

figures:
	$(PYTHON) tools/maintenance/run_experiment.py --experiment $(EXP) --profile full --figures $(if $(STUDY),--study $(STUDY),)

manuscript:
	mkdir -p manuscript/build
	cd manuscript && latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
	cd manuscript && latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=build supplementary.tex
	cp manuscript/build/main.pdf manuscript/manuscript.pdf
	cp manuscript/build/supplementary.pdf manuscript/supplementary-information.pdf

release:
	$(PYTHON) tools/maintenance/create_release_manifest.py --experiment $(EXP) --version $(VERSION)

clean:
	rm -rf manuscript/build/* artifacts/local/smoke .pytest_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

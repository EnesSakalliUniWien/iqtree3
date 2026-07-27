SHELL := /bin/bash
EXP ?= 003_invariant_mixture
PROFILE ?= smoke
STUDY ?=

.PHONY: help bootstrap doctor status map audit plot-check publication-figure-check figure-integration-check figure-artwork-check figure-equation-check figure-check figure-plan figure-render workflow-figure smoke check power-check rate-category-check data publication-figures experiment figures manuscript release clean

help:
	@echo "The sole supported workspace interface is ./satute"
	@echo "Run ./satute bootstrap once, then ./satute --help"

status:
	./satute status

bootstrap:
	./satute bootstrap

doctor:
	./satute doctor

map:
	./satute data --map

audit:
	./satute check audit

plot-check:
	./satute check publication-figures

publication-figure-check:
	./satute check publication-figures

figure-integration-check:
	./satute check figures

figure-artwork-check:
	./satute check figures

figure-equation-check:
	./satute check figures

figure-check: plot-check figure-integration-check figure-artwork-check figure-equation-check

figure-plan:
	./satute publication-figures --plan

figure-render:
	./satute publication-figures
	./satute workflow-figure

workflow-figure:
	./satute workflow-figure

smoke:
	./satute data --profile smoke

check:
	./satute check all

power-check:
	./satute check power

rate-category-check:
	./satute check rate-categories

data:
	./satute data --profile $(PROFILE)

publication-figures:
	./satute publication-figures

experiment:
	./satute experiment $(EXP) --profile $(PROFILE) $(if $(STUDY),--study $(STUDY),)

figures:
	./satute experiment $(EXP) --profile full --figures $(if $(STUDY),--study $(STUDY),)

manuscript:
	./satute manuscript

release:
	./satute release $(EXP) $(VERSION)

clean:
	./satute clean

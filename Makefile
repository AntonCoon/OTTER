PYTHON=python3

all:
	@mkdir -p bin
	@if test -d bin/snpEff; then echo snpEff exists; else echo Downloading snpEff... && cd bin \
		&& wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
		&& unzip snpEff_latest_core.zip; \
	fi
	@if test -d bin/snpEff/data/GRCh38.105; then echo reference genome found; else echo Downloading reference... && cd bin/snpEff/ \
		&& java -jar snpEff.jar download -v GRCh38.105; \
	fi
	@echo Seems all is fine...
	@echo
	@$(PYTHON) scripts/cli.py

test:
	# $(PYTHON) main.py --filter-hard example_data/BH_2/
	$(PYTHON) main.py --filter-medium example_data/BH_2/

.PHONY: clean
clean:
	rm -f tmp/*
	rm -rf example_data/BH_2/pipeline-out/
	rm -f bin

.PHONY: tags
tags:
	ctags -R . --languages=python,sh

cleanall: clean
	rm tags

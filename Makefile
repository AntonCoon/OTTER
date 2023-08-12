PYTHON=python3

all:
	$(PYTHON) scripts/cli.py

test:
	# $(PYTHON) main.py --filter-hard example_data/BH_2/
	$(PYTHON) main.py --filter-medium example_data/BH_2/

.PHONY: clean
clean:
	rm -f tmp/*
	rm -rf example_data/BH_2/pipeline-out/

.PHONY: tags
tags:
	ctags -R . --languages=python,sh

cleanall: clean
	rm tags

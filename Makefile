PYTHON=python3

all:
	$(PYTHON) scripts/cli.py

test:
	$(PYTHON) main.py --filte-hard example_data/BH_2/

.PHONY: clean
clean:
	rm -f tmp/*
	rm -f example_data/BH_2/pipeline-out/

.PHONY: tags
tags:
	ctags -R . --languages=python,sh

cleanall: clean
	rm tags

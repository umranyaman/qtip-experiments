PRED=python $(TS_HOME)/bin/predict.py
PRED_ARGS=

%.pred: %.ts
	$(PRED) --in $< --out $@ $(PRED_ARGS) --verbose --subsampling-fraction=0.3 --plot-all --subsampling-series 0.05,0.1,0.2,0.3,0.4,0.5,1.0 --subsampling-repl=5

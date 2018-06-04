
# -*- mode: Makefile -*-
#
# usage make -f THIS_FILE
#
DATATYPE := wals
LANGS := ../data/$(DATATYPE)/langs.json
FLIST := ../data/$(DATATYPE)/flist.json
OUTDIR := ../data/$(DATATYPE)

FNUM_FILE := $(OUTDIR)/fnum

SEED := --seed=2
PYTHON := nice -19 python

AUTO_OPTS := --norm_sigma=10.0 --gamma_scale=1.0 --burnin=500 --interval=5 --samples=500

LAST_MAKEFILE = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))


$(shell mkdir -p $(OUTDIR))

ifneq ("$(wildcard $(FNUM_FILE))","")
FEAT_NUM := $(shell cat $(FNUM_FILE))
else
FEAT_NUM := 0
endif
FEAT_MAX := $(shell expr $(FEAT_NUM) - 1)

$(FNUM_FILE) : $(FLIST)
	$(PYTHON) fnum.py $(FLIST) > $(FNUM_FILE)

# cv_auto_main FILE_PREFIX FID
define auto_main
$(1)_f$(2).json : $(OUTDIR)/langs.filled.json $(FLIST)
	$(PYTHON) run_autologistic.py $(SEED) $(AUTO_OPTS) --fid=$(2) $(OUTDIR)/langs.filled.json $(FLIST) $(1)_f$(2).json
AUTO_LIST += $(1)_f$(2).json
AUTO$(2)_LIST += $(1)_f$(2).json
endef

$(foreach f,$(shell seq 0 $(FEAT_MAX)), \
  $(eval $(call auto_main,$(OUTDIR)/auto,$(f))))


al_main : $(AUTO_LIST)

# recursive make; does not work?
al : $(FNUM_FILE)
	$(MAKE) -f $(LAST_MAKEFILE) al_main DATATYPE=$(DATATYPE) LANGS=$(LANGS) FLIST=$(FLIST) OUTDIR=$(OUTDIR) FNUM_FILE=$(FNUM_FILE) PYTHON="$(PYTHON)" SEED="$(SEED)" AUTO_OPTS="$(AUTO_OPTS)"



all : al

.phony : all al_main al

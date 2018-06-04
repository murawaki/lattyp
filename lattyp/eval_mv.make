# -*- mode: Makefile -*-
#
# usage make -f THIS_FILE CV=10
#
DATATYPE := wals
LANGS := ../data/$(DATATYPE)/langs.json
FLIST := ../data/$(DATATYPE)/flist.json
OUTDIR := ../data/$(DATATYPE)/cv

FNUM_FILE := ../data/$(DATATYPE)/fnum
CVMAP := $(OUTDIR)/langs.cvmap.json

NPB_KS := 50 100 250 500

CV := 10
CV_MAX := $(shell expr $(CV) - 1)

ifneq ("$(wildcard $(FNUM_FILE))","")
FEAT_NUM := $(shell cat $(FNUM_FILE))
else
FEAT_NUM := 0
endif
FEAT_MAX := $(shell expr $(FEAT_NUM) - 1)

MODEL_PREFIX := mda
MDA_KS := 50 100
TRAIN_OPTS := --maxanneal=100 --bias --norm_sigma=10.0 --gamma_scale=1.0 --resume_if
ITER := 500
SEED := --seed=2
PYTHON := nice -19 python
ITER_XZ := 100

AUTO_OPTS := --norm_sigma=10.0 --gamma_scale=1.0 --burnin=500 --interval=5 --samples=500

LAST_MAKEFILE = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))



$(shell mkdir -p $(OUTDIR))

$(FNUM_FILE) : $(FLIST)
	$(PYTHON) fnum.py $(FLIST) > $(FNUM_FILE)

$(CVMAP) : $(LANGS)
	$(PYTHON) -mmv.cv.make_cvmap $(SEED) $(LANGS) $(CVMAP) $(CV)


# cv_main FILE_PREFIX CV_IDX
define cv_main
$(1).cv$(2).json : $(LANGS) $(FLIST) $(CVMAP)
	$(PYTHON) -mmv.cv.hide $(LANGS) $(FLIST) $(1).cv$(2).json $(CVMAP) $(2)

HIDE_LIST += $(1).cv$(2).json

$(1).cv$(2).tsv : $(1).cv$(2).json $(FLIST)
	$(PYTHON) -mmv.json2tsv $(1).cv$(2).json $(FLIST) $(1).cv$(2).tsv

$(1).cv$(2).filled.tsv : $(1).cv$(2).tsv
	R --vanilla -f mv/impute_mca.r --args $(1).cv$(2).tsv $(1).cv$(2).filled.tsv

$(1).cv$(2).filled.json : $(1).cv$(2).filled.tsv $(FLIST)
	$(PYTHON) -mmv.tsv2json $(1).cv$(2).json $(1).cv$(2).filled.tsv $(FLIST) $(1).cv$(2).filled.json

FILLED_LIST += $(1).cv$(2).filled.json
MCR_LIST += $(1).cv$(2).filled.json
endef


$(foreach i,$(shell seq 0 $(CV_MAX)), \
  $(eval $(call cv_main,$(OUTDIR)/langs,$(i))))



##################################################
#                  BASIC                         #
##################################################

$(OUTDIR)/langs.random.eval : $(LANGS) $(FLIST)
	$(PYTHON) -mmv.cv.eval_mv $(SEED) --random $(LANGS) - $(FLIST) > $(OUTDIR)/langs.random.eval
EVALS_BASIC += $(OUTDIR)/langs.random.eval
EVALS += $(OUTDIR)/langs.random.eval

$(OUTDIR)/langs.freq.eval : $(LANGS) $(HIDE_LIST)
	sh -c 'for i in `seq 0 $(CV_MAX)`; do $(PYTHON) -mmv.cv.eval_mv $(SEED) --freq $(LANGS) $(OUTDIR)/langs.cv$${i}.json $(FLIST); done | perl -anle"\$$a+=\$$F[1];\$$b+=\$$F[2]; END{printf \"%f\\t%d\\t%d\\n\", \$$a / \$$b, \$$a, \$$b;}" > $(OUTDIR)/langs.freq.eval'

EVALS_BASIC += $(OUTDIR)/langs.freq.eval
EVALS += $(OUTDIR)/langs.freq.eval

$(OUTDIR)/langs.mcr.eval : $(LANGS) $(MCR_LIST)
	sh -c 'for i in `seq 0 $(CV_MAX)`; do $(PYTHON) -mmv.cv.eval_mv $(SEED) $(LANGS) $(OUTDIR)/langs.cv$${i}.filled.json -; done | perl -anle"\$$a+=\$$F[1];\$$b+=\$$F[2]; END{printf \"%f\\t%d\\t%d\\n\", \$$a / \$$b, \$$a, \$$b;}" > $(OUTDIR)/langs.mcr.eval'

EVALS_BASIC += $(OUTDIR)/langs.mcr.eval
EVALS += $(OUTDIR)/langs.mcr.eval



##################################################
#                   NPB                          #
##################################################

# cv_npb FILE_PREFIX CV_IDX K
define cv_npb_main
$(1).npb$(3).cv$(2).filled.json : $(1).cv$(2).tsv $(FLIST)
	mkdir -p $(1).npb$(3).cv$(2).filled.tsv.tmp \
	&& R --vanilla -f mv/impute_npb.r --args $(1).cv$(2).tsv $(1).npb$(3).cv$(2).filled.tmp $(3) \
	&& $(PYTHON) -mmv.tsv2json_merge $(1).cv$(2).json $(1).npb$(3).cv$(2).filled.tmp $(FLIST) $(1).npb$(3).cv$(2).filled.json \
	&& rm -r $(1).npb$(3).cv$(2).filled.tmp.* && rmdir $(1).npb$(3).cv$(2).filled.tsv.tmp

FILLED_LIST += $(1).npb$(3).cv$(2).filled.json
NPB_LIST += $(1).npb$(3).cv$(2).filled.json
NPB$(3)_LIST += $(1).npb$(3).cv$(2).filled.json
endef


$(foreach k,$(NPB_KS), \
   $(foreach i,$(shell seq 0 $(CV_MAX)), \
     $(eval $(call cv_npb_main,$(OUTDIR)/langs,$(i),$(k)))))


# cv_npb_eval FILE_PREFIX K
define cv_npb_eval
$(1).npb$(2).eval : $(NPB$(2)_LIST) $(LANGS)
	sh -c 'for i in `seq 0 $(CV_MAX)`; do $(PYTHON) -mmv.cv.eval_mv $(SEED) $(LANGS) $(1).npb$(2).cv$$$${i}.filled.json -; done | perl -anle"\$$$$a+=\$$$$F[1];\$$$$b+=\$$$$F[2]; END{printf \"%f\\t%d\\t%d\\n\", \$$$$a / \$$$$b, \$$$$a, \$$$$b;}" > $(1).npb$(2).eval'

EVALS += $(1).npb$(2).eval
EVALS_NPB += $(1).npb$(2).eval
endef

$(foreach k,$(NPB_KS), \
  $(eval $(call cv_npb_eval,$(OUTDIR)/langs,$(k))))



##################################################
#               AUTOLOGISTIC                     #
##################################################

# cv_auto_main FILE_PREFIX FID CV_IDX
define cv_auto_main
$(1)_f$(2)_CV$(3).json : $(OUTDIR)/langs.cv$(3).filled.json $(FLIST)
	$(PYTHON) run_autologistic.py $(SEED) $(AUTO_OPTS) --fid=$(2) --cv $(OUTDIR)/langs.cv$(3).filled.json $(FLIST) $(1)_f$(2)_CV$(3).json
AUTO_LIST += $(1)_f$(2)_CV$(3).json
AUTO$(2)_LIST += $(1)_f$(2)_CV$(3).json
endef

$(foreach f,$(shell seq 0 $(FEAT_MAX)), \
  $(foreach i,$(shell seq 0 $(CV_MAX)), \
    $(eval $(call cv_auto_main,$(OUTDIR)/auto,$(f),$(i)))))

al_main : $(AUTO_LIST)

# recursive make; does not work?
al : $(FNUM_FILE)
	$(MAKE) -f $(LAST_MAKEFILE) al_main DATATYPE=$(DATATYPE) LANGS=$(LANGS) FLIST=$(FLIST) OUTDIR=$(OUTDIR) CV=$(CV) FNUM_FILE=$(FNUM_FILE) PYTHON="$(PYTHON)" SEED="$(SEED)" AUTO_OPTS="$(AUTO_OPTS)"

$(OUTDIR)/auto.eval : $(AUTO_LIST)
	$(PYTHON) -mmv.merge_autologistic_results --cvout=$(OUTDIR)/auto.eval --detail=$(OUTDIR)/auto.results.json $(AUTO_LIST)

EVALS += $(OUTDIR)/auto.eval



##################################################
#                   MDA                          #
##################################################

# cv_mda_train FILE_PREFIX K CV_IDX
define cv_mda_train
$(1)_K$(2)_cv$(3).pkl : $(OUTDIR)/langs.cv$(3).filled.json $(FLIST)
	$(PYTHON) train_mda.py $(SEED) --iter=$(ITER) --cv --K=$(2) --output=$(1)_K$(2)_cv$(3).pkl $(OUTDIR)/langs.cv$(3).filled.json $(TRAIN_OPTS) $(FLIST) >> $(1)_K$(2)_cv$(3).log 2>&1 && cp $(1)_K$(2)_cv$(3).pkl.final $(1)_K$(2)_cv$(3).pkl

MDAS += $(1)_K$(2)_cv$(3).pkl
endef

$(foreach k,$(MDA_KS), \
  $(foreach i,$(shell seq 0 $(CV_MAX)), \
    $(eval $(call cv_mda_train,$(OUTDIR)/$(MODEL_PREFIX),$(k),$(i)))))


# cv_mda_sample_auto FILE_PREFIX K CV_IDX
define cv_mda_sample_auto
$(1)_K$(2)_cv$(3).xz.json.bz2 : $(1)_K$(2)_cv$(3).pkl
	$(PYTHON) sample_auto.py $(SEED) --a_repeat=5 --iter=$(ITER_XZ) $(1)_K$(2)_cv$(3).pkl - | bzip2 -c > $(1)_K$(2)_cv$(3).xz.json.bz2

SAMPLE_AUTO += $(1)_K$(2)_cv$(3).xz.json.bz2
SAMPLE_XZ_$(2) += $(1)_K$(2)_cv$(3).xz.json.bz2
SAMPLE_XZ_$(2)_$(3) += $(1)_K$(2)_cv$(3).xz.json.bz2

$(1)_K$(2)_cv$(3).xz.merged.json : $(1)_K$(2)_cv$(3).xz.json.bz2
	$(PYTHON) convert_auto_xz.py --burnin=0 --update --input=$(1)_K$(2)_cv$(3).xz.json.bz2 $(OUTDIR)/langs.cv$(3).filled.json $(FLIST) > $(1)_K$(2)_cv$(3).xz.merged.json

SAMPLE_XZ_MERGED += $(1)_K$(2)_cv$(3).xz.merged.json
SAMPLE_XZ_MERGED_$(2) += $(1)_K$(2)_cv$(3).xz.merged.json
endef

# cv_mda_eval FILE_PREFIX K
define cv_mda_eval
$(1)_K$(2).eval : $(SAMPLE_XZ_MERGED_$(2))
	sh -c 'for i in `seq 0 $(CV_MAX)`; do $(PYTHON) -mmv.cv.eval_mv $(SEED) $(LANGS) $(1)_K$(2)_cv$$$${i}.xz.merged.json -; done | perl -anle"\$$$$a+=\$$$$F[1];\$$$$b+=\$$$$F[2]; END{printf \"%f\\t%d\\t%d\\n\", \$$$$a / \$$$$b, \$$$$a, \$$$$b;}" > $(1)_K$(2).eval'

EVALS_MDA += $(1)_K$(2).eval
EVALS += $(1)_K$(2).eval
endef

$(foreach k,$(MDA_KS), \
  $(foreach i,$(shell seq 0 $(CV_MAX)), \
     $(eval $(call cv_mda_sample_auto,$(OUTDIR)/$(MODEL_PREFIX),$(k),$(i)))))
$(foreach k,$(MDA_KS), \
  $(eval $(call cv_mda_eval,$(OUTDIR)/$(MODEL_PREFIX),$(k))))



all : $(EVALS)

hidden : $(HIDE_LIST)

filled : $(FILLED_LIST)

mcr : $(OUTDIR)/langs.mcr.eval

basic : $(EVALS_BASIC)

npb : $(EVALS_NPB)

auto : $(OUTDIR)/auto.eval

mda_model : $(MDAS)

mda : $(EVALS_MDA)

clean :
	rm -f $(OUTDIR)/langs*
	rmdir --ignore-fail-on-non-empty $(OUTDIR)

.PHONY : all clean hidden filled mcr basic npb al_main al auto mda_model mda

# Builds the post-processing programs in tools/ and runs the regression checks.
#
#   make              compile every tools/*.cpp into build/bin/<name>
#   make pp_meas_ACs  compile one tool (any name from tools/)
#   make test         frozen-file and hygiene checks, then run every tool on the
#                     stored reference dumps and compare with a tolerance
#   make test-exact   the same, but every output must match byte for byte
#   make clean        remove build/ (and nothing else)
#   make help         this list
#
# Variables you can set on the command line:
#
#   CXX=...           compiler (default g++, as in v2.0.0)
#   LEGACY_CXXFLAGS=  flags for the tools; empty by default, which is the
#                     plain "g++ -o name name.cpp" used for the published
#                     results and for the stored reference answers
#   BUILD=dir         where binaries and regression runs go (default build)
#   REGRESS_SKIP=...  comma-separated tools left out of the regression run.
#                     Default pp_meas_Contours: it needs 4 to 5.5 GB of memory
#                     and writes two files of 300 to 400 MB per reference set.
#                     REGRESS_SKIP= runs it too.
#   RTOL=, ATOL=      tolerances for "make test" (default 1e-6 and 0). Keep
#                     ATOL at 0: the outputs are in SI units and some columns
#                     are of order 1e-13 (s) or 1e-28 (J), so any fixed
#                     absolute tolerance would hide changes there.
#
# The helper targets of the old Makefile (collect, grep, update, q) are in
# scripts/rundir.mk and run inside a simulation directory.

# Default compiler: g++, like the manual build of v2.0.0. make's own default
# is c++ on some systems, so replace it unless CXX was given by the user.
ifeq ($(origin CXX),default)
CXX := g++
endif

# Keep empty to reproduce the stored reference outputs byte for byte.
LEGACY_CXXFLAGS :=

BUILD        ?= build
BINDIR       := $(BUILD)/bin
REGRESS_WORK := $(BUILD)/regress
REGRESS_SKIP ?= pp_meas_Contours
RTOL         ?= 1e-6
ATOL         ?= 0
PYTHON       ?= python3

LEGACY_SRC := $(wildcard tools/*.cpp)
LEGACY     := $(basename $(notdir $(LEGACY_SRC)))
LEGACY_BIN := $(addprefix $(BINDIR)/,$(LEGACY))

# sha256 check of the stored reference files (GNU or BSD tools).
SHA256_CHECK := $(shell command -v sha256sum >/dev/null 2>&1 && echo "sha256sum -c --quiet" || echo "shasum -a 256 -c --quiet")

# Reference sets with a "mini" regression directory.
REGRESS_CASES := $(patsubst %/DUMP_CONTENT.sha256,%,$(wildcard \
    tests/reference/R3_channel_mini/mini/DUMP_CONTENT.sha256 \
    tests/reference/R2_examples_mini/mini/DUMP_CONTENT.sha256))

.PHONY: all test test-exact checks clean help $(LEGACY)

all: $(LEGACY_BIN)

$(BINDIR)/%: tools/%.cpp
	@mkdir -p $(BINDIR)
	$(CXX) $(LEGACY_CXXFLAGS) -o $@ $<

# "make pp_meas_ACs" builds that one tool.
$(LEGACY): %: $(BINDIR)/%

# Without --exact, pp_meas_wallTemp is not run on a set without a wall dump:
# it would work on values that were never set (see run_legacy_tools.sh).
test: all checks
	@$(MAKE) --no-print-directory regress SKIP_MISSING_INPUT=1 \
	    COMPARE_OPTS="--rtol $(RTOL) --atol $(ATOL)"

test-exact: all checks
	@$(MAKE) --no-print-directory regress SKIP_MISSING_INPUT=0 COMPARE_OPTS="--exact"

checks:
	$(PYTHON) tests/regress/check_frozen.py
	$(PYTHON) tests/regress/check_hygiene.py

# Runs each tool on every reference set and compares the outputs; checks
# all sets before reporting a failure.
.PHONY: regress
regress:
	@if [ -z "$(strip $(REGRESS_CASES))" ]; then \
	    echo "no reference sets found under tests/reference/*/mini"; exit 1; \
	fi
	@status=0; n=0; \
	for ref in $(REGRESS_CASES); do \
	    name=$$(basename $$(dirname $$ref)); \
	    echo "== $$name"; n=$$((n+1)); \
	    if ! (cd $$ref && $(SHA256_CHECK) SHA256SUMS); then \
	        echo "stored reference files differ from $$ref/SHA256SUMS"; status=1; continue; \
	    fi; \
	    if SKIP_TOOLS="$(REGRESS_SKIP)" SKIP_MISSING_INPUT=$(SKIP_MISSING_INPUT) \
	            bash tests/regress/run_legacy_tools.sh \
	            $$ref $(BINDIR) $(REGRESS_WORK)/$$name; then \
	        $(PYTHON) tests/regress/compare_outputs.py $$ref/outputs \
	            $(REGRESS_WORK)/$$name/tools --skip "$(REGRESS_SKIP)" \
	            $(COMPARE_OPTS) || status=1; \
	    else \
	        status=1; \
	    fi; \
	done; \
	echo "checked $$n reference sets"; \
	exit $$status

clean:
	@case "$(BUILD)" in ""|.|./|..|../|/) echo "not removing BUILD=$(BUILD)"; exit 1 ;; esac
	rm -rf $(BUILD)

# Prints the comment block at the top of this file.
help:
	@sed -n '/^$$/q; s/^# \{0,1\}//p' $(firstword $(MAKEFILE_LIST))

# Helper targets for a simulation directory, kept from the Makefile of
# v2.0.0. Run them from inside the directory that holds the dump files:
#
#   cd my_run
#   make -f /path/to/ScatteringKernels/scripts/rundir.mk update
#
#   collect  move the tool outputs (*.txt *.dat *.trj *.xyz, but not
#            Specification.dat or data.dat) into DataDir/
#   grep     print the number of frames in dump_meas_gas.lammpstrj
#   update   write that number into line 1 of Specification.dat, which
#            the tools read as the number of frames to process
#   q        list your jobs in the Slurm queue
#
# Change the names with DATADIR=..., DUMP=... or SPEC=... on the command line.

DATADIR ?= DataDir
DUMP    ?= dump_meas_gas.lammpstrj
SPEC    ?= Specification.dat

.PHONY: collect grep update q

collect:
	@mkdir -p $(DATADIR)
	@for f in *.txt *.dat *.trj *.xyz; do \
	    case "$$f" in $(SPEC)|data.dat) continue ;; esac; \
	    if [ -f "$$f" ]; then mv -f "$$f" $(DATADIR)/ && echo "moved $$f"; fi; \
	done

grep:
	@grep -c 'ITEM: TIMESTEP' $(DUMP)

# Writes a temporary file and moves it back, so it works with BSD and GNU sed.
update:
	@n=`grep -c 'ITEM: TIMESTEP' $(DUMP)` && \
	echo "nTimeSteps = $$n" && \
	sed "1s/.*/$$n                      # nTimeSteps/" $(SPEC) > $(SPEC).tmp && \
	mv $(SPEC).tmp $(SPEC)

q:
	squeue -u $(USER)

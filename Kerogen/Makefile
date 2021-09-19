# 1. Make attemps to build the first target it finds in the Makefile, the default target
# 2. The data files are a product not only of our text files (e.g. dump files) but the script, 
# that processes the text files and creates the data files.
# 3. Updating a subset of the files in the pipeline triggers rerunning the appropriate downstream steps.
# 4. Use the wildcard "%" as a placeholder in targets and dependencies, 
# while use the special variable "$*" to refer to matching sets of files in actions.
# 5. Define variables by assigning values to names, reference variable using $(...)
# 6. Make is actually a small programming language with many built-in functions. 
# 7. Use wildcard "*" function to get lists of files matching a pattern. 
# 8. Use patsubst "%" function to rewrite file names. 
# 9. Document Makefiles by adding specially-formatted comments and a target to extract and format them. 

# $@ means "the target of the current rule"
# $^ means "all the dependencies of the current rule"
# $< refers to "the first dependency of the current rule"

HOME    =/work/d114/d114/yichongc/Kerogen/Kerogen_On
INCL    =-I${HOME}/incl
SRC     =${HOME}/src
BIN     =${HOME}/bin
DataDir =${HOME}/DataDir

.SUFFIXES: .cpp
# It informs make that you will be using this special suffixes to make your own rules.
# Implicit rule for the construction of .o (object) files out of .cpp (source files)
%.o : ${SRC}/%.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ $<
# The -c flag says to generate the objective file. 
# The -o $@ says to put the output of the compilation in the file named on the left side of the :	


# Start to compile the codes
flowRate : ${SRC}/pp_meas_flowRate.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

profiles : ${SRC}/pp_meas_bin.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

ACs : ${SRC}/pp_meas_ACs.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

EAC : ${SRC}/pp_meas_EAC.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

FVDF : ${SRC}/pp_meas_FVDF.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

PVDF : ${SRC}/pp_meas_PVDF.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

PIVDF : ${SRC}/pp_meas_PIVDF.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

FAngular :${SRC}/pp_meas_FAngular.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

PAngular :${SRC}/pp_meas_PAngular.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

MagV : ${SRC}/pp_meas_MagV.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Time : ${SRC}/pp_meas_time.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

KerogenTime : ${SRC}/pp_kerogen_time.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Space : ${SRC}/pp_meas_space.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Traj : ${SRC}/pp_meas_traj.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Collisions : ${SRC}/pp_meas_collisions.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

KerogenEquil : ${SRC}/pp_equil_kerogen.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

KerogenMeas : ${SRC}/pp_meas_kerogen.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Correlation : ${SRC}/pp_meas_correlation.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

# ------------------------------------------------------------
# Phony target, that it does not build anything

## PoiseData	: Compile the source codes for Poiseuille Flow simulation.

.PHONY : AccommodationCoefficients
AccommodationCoefficients : ACs EAC

.PHONY : KerogenData
KerogenData :  FVDF PVDF ACs EAC profiles Collisions FAngular PAngular Traj Time Space KerogenMeas flowRate Correlation 
	
.PHONY : PlatinumData
PlatinumData : FVDF PVDF ACs EAC profiles Collisions FAngular PAngular Traj Time Space


## collect	: Collect the auto-generated files
.PHONY : collect
collect :
	mv -f *.txt ${DataDir}
	mv -f *.dat ${DataDir}
	mv ${DataDir}/nTimeSteps.dat ./
	mv ${DataDir}/data.dat ./
	mv ${DataDir}/Height.dat ./

## clean	: Remove the auto-generated files
.PHONY : clean
clean :
	rm -f *.txt
	rm -f *.dat
	rm -f *.o

## q		: Monitoring jobs
.PHONY : q
q : 
	squeue -u $(USER)

## grep		: Extract the timesteps in the dumpfile
.PHONY : grep
grep : 
	grep -o 'TIMESTEP' dump_meas_gas.lammpstrj | wc -l

# sed stands for 'stream editor'. sed reads in some text, does some filtering, and writes out the filtered text. 
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<


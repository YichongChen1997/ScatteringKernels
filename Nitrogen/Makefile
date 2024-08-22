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

HOME    =/work/d114/d114/yichongc/Nitrogen_300K/Nitrogen_0.00MPa
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
FlowRate : ${SRC}/pp_meas_FlowRate.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

ACs : ${SRC}/pp_meas_ACs.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Bins : ${SRC}/pp_meas_Bins.o
	$(CXX) -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

FVDF : ${SRC}/pp_VDF_F.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

PrVDF : ${SRC}/pp_VDF_Pr.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

PiVDF : ${SRC}/pp_VDF_Pi.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

FAngular :${SRC}/pp_Angular_F.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

PAngular :${SRC}/pp_Angular_P.o
	${CXX} -o $@ $^
	rm -r $^
	./$@
	mv $@ ${BIN}

ResidenceTime : ${SRC}/pp_meas_ResidenceTime.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Displacement : ${SRC}/pp_meas_Displacement.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Traj : ${SRC}/pp_meas_Traj.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

Collisions : ${SRC}/pp_meas_Collisions.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

VirtualLayer : ${SRC}/pp_meas_VirtualLayer.o
	${CXX} -o $@ $^
	rm -f $^
	./$@
	mv $@ ${BIN}

# ------------------------------------------------------------
# Phony target, that it does not build anything

## PoiseData	: Compile the source codes for Poiseuille Flow simulation.
.PHONY : PoiseData
PoiseData : FlowRate Bins

.PHONY : AccommodationCoefficients
AccommodationCoefficients : ACs Collisions

.PHONY : ScatteringData
ScatteringData : FAngular PAngular ACs Bins Collisions Displacement FlowRate ResidenceTime Traj VirtualLayer FVDF PiVDF PrVDF

## collect	: Collect the auto-generated files
.PHONY : collect
collect :
	mv -f *.txt ${DataDir}
	mv -f *.dat ${DataDir}
	mv -f *.trj ${DataDir}
	mv -f *.xyz ${DataDir}
	mv ${DataDir}/Specification.dat ./
	mv ${DataDir}/data.dat ./

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


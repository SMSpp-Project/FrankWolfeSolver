##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of FrankWolfeSolver                                             #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by FrankWolfeSolver, i.e., core SMS++ only: the         #
#   :Solver acting as Linear Minimization Oracles for the sub-Block (and     #
#   the :Solver used for comparison) are obtained through the Solver         #
#   factory and the CDASolver interface, hence they are not compile-time     #
#   dependencies (exactly as for LagrangianDualSolver and its inner          #
#   Solver).                                                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by FrankWolfeSolver.                                                 #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(FWSlvSDR)    = the directory where the source is               #
#                                                                            #
#   Output: $(FWSlvOBJ)    = the final object(s) / library                   #
#           $(FWSlvH)      = the .h files to include                         #
#           $(FWSlvINC)    = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

FWSlvOBJ = $(FWSlvSDR)/obj/FrankWolfeSolver.o

FWSlvINC = -I$(FWSlvSDR)/include

FWSlvH   = $(FWSlvSDR)/include/FrankWolfeSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(FWSlvOBJ) $(FWSlvSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(FWSlvSDR)/obj/FrankWolfeSolver.o: \
	$(FWSlvSDR)/src/FrankWolfeSolver.cpp \
	$(FWSlvSDR)/include/FrankWolfeSolver.h $(SMS++OBJ)
	$(CC) -c $(FWSlvSDR)/src/FrankWolfeSolver.cpp -o $@ \
	$(FWSlvINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################

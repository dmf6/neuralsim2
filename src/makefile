SHELL = /bin/sh

CC = g++
CFLAGS= 
LIBS = 
CPPFLAGS = -g -O2 -Wall

SUN_DISTRO = /usr/local/bin/sundials-config

VPATH=%.h ./include
VPATH=%.o ./obj

OBJDIR = ./obj
INCLUDE_DIR = ./include
INCLUDES  := $(addprefix -I,$(INCLUDE_DIR))

objects = $(addprefix $(OBJDIR)/, electrode.o mm.o neuron.o neuronmodel.o rk.o synapse.o simplecell.o random.o abneuron.o gap.o CVodeSolver.o chemsyn.o rk65n.o)

MY_APPS = neuralsim

#determine the preprocessor and linker flags required to link the
#SUNDIALS libraries to the user application by invoking, from within
#the Makefile, the executable script sundials-config

all:
	sun_cppflags=`eval "${SUN_DISTRO} -m cvode -t s -l c -s cppflags"`; \
	sun_ldflags=`eval "${SUN_DISTRO} -m cvode -t s -l c -s libs"`;      \
	for i in ${MY_APPS} ; do                                                                \
           echo "--- Making $${i} ---" ;                                                         \
           eval "make SUN_CPPFLAGS='$${sun_cppflags}' SUN_LDFLAGS='$${sun_ldflags}' $${i}";      \
	done

${MY_APPS} : $(objects)
	${CC} -o ${MY_APPS} ${objects} ${SUN_LDFLAGS} ${LIBS}

$(OBJDIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) ${SUN_CPPFLAGS}  ${INCLUDES} $< -o $@


$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	rm -f ${MY_APPS}
	rm -f ${objects}

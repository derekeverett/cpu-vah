DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic
DIR_H          = $(DIR_MAIN)include/
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O5
FLOWTRACE =
OPTIONS =
LINK_OPTIONS = -link -L/home/everett.165/libconfig-1.5/lib/.libs -lconfig -L/home/everett.165/googletest-master/googletest/mybuild/ -lgtest
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = nvcc
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest
INCLUDES = -I rhic/rhic-core/src/include -I rhic/rhic-harness/src/main/include -I rhic/rhic-trunk/src/include -I rhic/rhic-harness/src/include -I /home/everett.165/libconfig-1.5/lib/ -I /home/everett.165/googletest-master/googletest/include/ -I freezeout

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)

EXE =\
cpu-vah

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic/rhic-core rhic/rhic-harness rhic/rhic-trunk -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:

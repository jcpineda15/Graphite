CC		= g++
CFLAGS	= -std=c++14 -O2 -Wall
LDFLAGS = 
LIBS	= -lz -lpthread
PROGRAM = ./graphite

SOURCE_DIR = src
HEADER_DIR = include
OUT_DIR = .
OBJ_DIR = obj
DEPEND_DIR = depend

SOURCES = $(wildcard $(SOURCE_DIR)/*.cpp)
HEADERS = $(wildcard $(HEADER_DIR)/*.h)

SOURCE_NAMES = $(notdir $(SOURCES))
OBJS = $(addprefix $(OBJ_DIR)/,$(SOURCE_NAMES:.cpp=.o))
DEPENDS = $(addprefix $(DEPEND_DIR)/,$(SOURCE_NAMES:.cpp=.depend))

.PHONY: all
all: $(DEPEND_DIR) $(OBJ_DIR) $(DEPENDS) $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $(PROGRAM) $(LIBS)

$(DEPENDS): | $(DEPEND_DIR)

$(DEPEND_DIR)/%.depend: $(SOURCE_DIR)/%.cpp $(HEADERS)
	@echo "generating $@"
	@$(CC) $(CFLAGS) $(LIBS) -MM $< > $@

$(DEPEND_DIR):
	mkdir -p -v $(DEPEND_DIR)

$(OBJS): | $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) $(LIBS) -c $^ -o $@

$(OBJ_DIR):
	mkdir -p -v $(OBJ_DIR)

ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEPENDS)
endif

.PHONY : clean
clean:
	rm -r $(OBJ_DIR)/*.o $(PROGRAM)

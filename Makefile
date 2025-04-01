DIR_SRC := src
DIR_INC := include
DIR_BIN := .build
DIR_LIB := lib

DIR_GUARD = mkdir -p $(@D)

CC := clang
LIBS := -lm -lsndfile
PKG_CONFIG := #libpipewire-0.3
CFLAGS := -O3 -Wextra -Wpedantic -std=gnu23 -I$(DIR_INC)
DEFNS := -DDSP_SAMPLE_FLOAT

TARGET := $(DIR_LIB)/libdsp.so
HEADERS := $(wildcard $(DIR_INC)/dsp/*.h)
BINARIES := $(patsubst $(DIR_SRC)/%.c, $(DIR_BIN)/%.o, $(wildcard $(DIR_SRC)/*.c))

$(TARGET): $(BINARIES)
	$(DIR_GUARD)
	$(CC) -fPIC $^ -shared $(shell pkg-config --libs $(PKG_CONFIG)) $(LIBS) -o $@

$(DIR_BIN)/%.o: $(DIR_SRC)/%.c $(HEADERS)
	$(DIR_GUARD)
	$(CC) $(DEFNS) $(CFLAGS) $(shell pkg-config --cflags $(PKG_CONFIG)) -fPIC -c $< -o $@ 

$(DIR_BIN)/test: test/test.c $(TARGET)
	$(CC) $(DEFNS) $(CFLAGS) $< $(shell pkg-config --libs $(PKG_CONFIG)) $(LIBS) -L$(DIR_LIB) -ldsp -o $@

.PHONY: test test-% clean all

test: $(DIR_BIN)/test
	LD_LIBRARY_PATH=$(DIR_LIB):${LD_LIBRARY_PATH} $^

test-%: $(DIR_BIN)/test
	LD_LIBRARY_PATH=$(DIR_LIB):${LD_LIBRARY_PATH} $^ $*

clean: 
	rm -rf $(DIR_BIN)
	rm -rf $(DIR_LIB)

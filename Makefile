LIB_NAME := dsp

DIR_SRC := src
DIR_INC := include
DIR_BIN := .build
DIR_LIB := lib
DIR_IST := ${HOME}/.local

DIR_GUARD = mkdir -p $(@D)

CC := clang
LIBS := -lm -lsndfile
PKG_CONFIG := #libpipewire-0.3
CFLAGS := -O3 -Wextra -Wpedantic -std=gnu23 -I$(DIR_INC)
DEFNS := -DDSP_SAMPLE_FLOAT

TARGET := $(DIR_LIB)/lib$(LIB_NAME).so
HEADERS := $(wildcard $(DIR_INC)/$(LIB_NAME)/*.h)
BINARIES := $(patsubst $(DIR_SRC)/%.c, $(DIR_BIN)/%.o, $(wildcard $(DIR_SRC)/*.c))

$(TARGET): $(BINARIES)
	$(DIR_GUARD)
	$(CC) -fPIC $^ -shared $(shell pkg-config --libs $(PKG_CONFIG)) $(LIBS) -o $@

$(DIR_BIN)/%.o: $(DIR_SRC)/%.c $(HEADERS)
	$(DIR_GUARD)
	$(CC) $(DEFNS) $(CFLAGS) $(shell pkg-config --cflags $(PKG_CONFIG)) -fPIC -c $< -o $@ 

$(DIR_BIN)/test: test/test.c $(TARGET)
	$(CC) $(DEFNS) $(CFLAGS) $< $(shell pkg-config --libs $(PKG_CONFIG)) $(LIBS) -L$(DIR_LIB) -l$(LIB_NAME) -o $@

.PHONY: test test-% clean install uninstall

install: $(TARGET)
	mkdir -p $(DIR_IST)/$(DIR_LIB)
	cp $(TARGET) $(DIR_IST)/$(DIR_LIB)/lib$(LIB_NAME).so
	mkdir -p $(DIR_IST)/$(DIR_INC)
	cp -r $(DIR_INC)/$(LIB_NAME) $(DIR_IST)/$(DIR_INC)

uninstall: 
	rm $(DIR_IST)/$(DIR_LIB)/lib$(LIB_NAME).so
	rm -r $(DIR_IST)/$(DIR_INC)/$(LIB_NAME)

test: $(DIR_BIN)/test
	LD_LIBRARY_PATH=$(DIR_LIB):${LD_LIBRARY_PATH} $^

test-%: $(DIR_BIN)/test
	LD_LIBRARY_PATH=$(DIR_LIB):${LD_LIBRARY_PATH} $^ $*

clean: 
	rm -rf $(DIR_BIN)
	rm -rf $(DIR_LIB)

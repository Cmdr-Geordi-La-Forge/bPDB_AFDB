CC       = gcc
CFLAGS   = -O3 -Wall -Wextra -std=c99 -march=native -Iinclude -Isrc
DEPFLAGS = -MMD -MP

BUILD_DIR = build
EXEC      = $(BUILD_DIR)/test_encoder

SRCS = \
	src/encoder.c \
	src/decoder.c \
	src/binary/meta_codec.c \
	src/binary/coord_codec.c \
	src/chem/dict.c \
	src/chem/struct_conf.c \
	src/cif/cif_lexer.c \
	src/cif/cif_parser.c \
	src/cif/cif_writer.c \
	src/cif/afdb_strings.c \
	src/util/str_util.c \
	tests/test_roundtrip.c

OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRCS))
DEPS = $(OBJS:.o=.d)

.PHONY: all test clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

test: $(EXEC)
	./$(EXEC) tests/fixtures/AF-A0A011QK89-F1.cif

clean:
	rm -rf $(BUILD_DIR) output.bpdb decoded_AF-*.cif

-include $(DEPS)

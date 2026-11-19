CC       = gcc
CFLAGS   = -O3 -Wall -Wextra -std=c99 -march=native -Iinclude -Isrc
DEPFLAGS = -MMD -MP

BUILD_DIR = build
LIB_A     = $(BUILD_DIR)/libbpdb.a
EXEC      = $(BUILD_DIR)/test_encoder
HARNESS   = $(BUILD_DIR)/bpdb_harness

LIB_SRCS = \
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
	src/util/str_util.c

HARNESS_SRCS = \
	tools/harness.c \
	tools/tar_stream.c \
	tools/gz_inflate.c \
	tools/diff_bucket.c

LIB_OBJS     = $(patsubst %.c,$(BUILD_DIR)/%.o,$(LIB_SRCS))
TEST_OBJS    = $(BUILD_DIR)/tests/test_roundtrip.o
HARNESS_OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(HARNESS_SRCS))
DEPS         = $(LIB_OBJS:.o=.d) $(TEST_OBJS:.o=.d) $(HARNESS_OBJS:.o=.d)

HARNESS_CFLAGS  = -Itools
HARNESS_LDFLAGS = -lpthread
ifeq ($(USE_LIBDEFLATE),1)
	HARNESS_CFLAGS  += -DBPDB_USE_LIBDEFLATE
	HARNESS_LDFLAGS += -ldeflate
else
	HARNESS_LDFLAGS += -lz
endif

.PHONY: all test harness clean

all: $(EXEC) $(HARNESS)

$(LIB_A): $(LIB_OBJS)
	ar rcs $@ $^

$(EXEC): $(TEST_OBJS) $(LIB_A)
	$(CC) $(CFLAGS) -o $@ $^

$(HARNESS): $(HARNESS_OBJS) $(LIB_A)
	$(CC) $(CFLAGS) -o $@ $^ $(HARNESS_LDFLAGS)

$(BUILD_DIR)/tools/%.o: tools/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(HARNESS_CFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

test: $(EXEC)
	./$(EXEC) tests/fixtures/AF-A0A011QK89-F1.cif

harness: $(HARNESS)

clean:
	rm -rf $(BUILD_DIR) output.bpdb decoded_AF-*.cif

-include $(DEPS)

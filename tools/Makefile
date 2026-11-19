INGEST    = $(BUILD_DIR)/bpdb_ingest
INDEX     = $(BUILD_DIR)/bpdb_index

INGEST_SRCS = \
	tools/bpdb_ingest.c \
	tools/pack_writer.c \
	tools/gcs_stream.c \
	tools/gz_inflate.c

INDEX_SRCS = \
	tools/bpdb_index.c

INGEST_OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(INGEST_SRCS))
INDEX_OBJS  = $(patsubst %.c,$(BUILD_DIR)/%.o,$(INDEX_SRCS))
DEPS        = $(LIB_OBJS:.o=.d) $(TEST_OBJS:.o=.d) $(HARNESS_OBJS:.o=.d) $(INGEST_OBJS:.o=.d) $(INDEX_OBJS:.o=.d)

INGEST_LDFLAGS = $(HARNESS_LDFLAGS) -lcurl -larchive

.PHONY: all test harness ingest index clean

all: $(EXEC) $(HARNESS) $(INGEST) $(INDEX)

$(INGEST): $(INGEST_OBJS) $(LIB_A)
	$(CC) $(CFLAGS) -o $@ $^ $(INGEST_LDFLAGS)

$(INDEX): $(INDEX_OBJS) $(LIB_A)
	$(CC) $(CFLAGS) -o $@ $^

ingest: $(INGEST)
index: $(INDEX)

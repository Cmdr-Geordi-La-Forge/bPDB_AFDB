CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -march=native

SRCS = test.c encoder.c decoder.c dict.c utils.c pdb_utils.c afdb_strings.c
OBJS = $(SRCS:.c=.o)
EXEC = test_encoder

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c bpdb.h afdb_strings.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(EXEC)
	./$(EXEC)

clean:
	rm -f $(OBJS) $(EXEC) output.bpdb decoded_AF-*.cif

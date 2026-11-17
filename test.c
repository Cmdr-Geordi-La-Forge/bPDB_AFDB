#include "bpdb.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    const char* filename = "AF-A0A011QK89-F1.cif";
    FILE* f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Error: Could not open '%s' for testing.\n", filename);
        return 1;
    }

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char* cif_data = malloc(fsize + 1);
    if (!cif_data) {
        fclose(f);
        return 1;
    }
    
    if (fread(cif_data, 1, fsize, f) != (size_t)fsize) {
        fprintf(stderr, "Error: Failed to read file data.\n");
        free(cif_data); fclose(f); return 1;
    }
    fclose(f);
    cif_data[fsize] = '\0';

    uint8_t* out_buf = NULL;
    size_t bin_size = encode_afdb_to_binary(cif_data, fsize, &out_buf);

    if (bin_size > 0 && out_buf) {
        printf("Success: Encoded %ld byte mmCIF to %zu byte binary stream.\n", fsize, bin_size);
        FILE* out = fopen("output.bpdb", "wb");
        if (out) {
            fwrite(out_buf, 1, bin_size, out);
            fclose(out);
            printf("Saved encoded binary to 'output.bpdb'\n");
            
            FILE* dec_out = fopen("decoded_AF-A0A011QK89-F1.cif", "w");
            if (dec_out) {
                decode_binary_to_afdb(out_buf, bin_size, dec_out);
                fclose(dec_out);
                printf("Successfully decoded back to 'decoded_AF-A0A011QK89-F1.cif'\n");
            }
        }
    } else {
        fprintf(stderr, "Failure: Encoding returned 0 bytes.\n");
    }

    free(cif_data);
    free(out_buf);
    return 0;
}

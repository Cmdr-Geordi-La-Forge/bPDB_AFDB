# Binary PDB (bPDB) v2 Format Specification

**bPDB v2** is a lossless, two-tier binary container format for AlphaFold Database (AFDB) and ModelCIF protein structures. It achieves ~24× lossless compression over uncompressed mmCIF while supporting zero-allocation, `mmap`-friendly SIMD geometry decoding (`x86-64 SSE4.1` and `AArch64 NEON`) and exact byte-for-byte mmCIF roundtrip reconstruction.

> **Deprecation Notice:** The legacy v0 PDB-centric format (`version == 0`, 8-byte aligned chain clusters, 13/15-bit sequential deltas, CONECT stream) is deprecated and superseded by v2 (`version == 2`).

---

## 1. Two-Tier File Architecture

Every bPDB v2 stream is **little-endian** and partitioned into a 4-byte-aligned **Hot Geometry Tier** followed immediately by a tightly packed **Cold Metadata Tier**:

```text
+----------------------------------------------------------------------+ offset 0
| Fixed Global Header (BpdbV2Header)                         64 bytes  |
+----------------------------------------------------------------------+ offset 64
| HOT GEOMETRY TIER (4-byte aligned, zero-copy mmap, SIMD-ready)       |
|   1. Sequence Array        ASCII 'A'..'Z'         align4(total_res)  |
|   2. pLDDT Array           uint16_t (scaled x100) align4(total_res*2)|
|   3. 12-Bit Delta Stream   9 bytes per 2 atoms    coord_padded_size()|
|   4. Escape Atom Table     EscapeAtom[esc_count]  escape_count * 12  |
+----------------------------------------------------------------------+ offset = cold_offset
| COLD METADATA TIER (lossless mmCIF reconstruction, tightly packed)   |
|   1. _atom_site Widths     uint8_t[8]             8 bytes            |
|   2. Flag-Gated Sections   Compact varlen / TLV / preset payloads    |
+----------------------------------------------------------------------+ EOF
```

* **Hot Geometry Tier (`[0 .. cold_offset)`)**: Self-contained for structural biology, geometric analysis, and ML pipelines (UniProt accession, fragment index, 1-letter sequence, per-residue pLDDT, and 3D heavy-atom coordinates). Opened with zero heap allocations via `bpdb_open_hot()` and decoded directly into `float4[x, y, z, 0]` via `bpdb_decode_coords_f32()`.
* **Cold Metadata Tier (`[cold_offset .. EOF)`)**: Stores formatting widths, global pLDDT, entity description, delta-encoded secondary structure (`_struct_conf`), default-elided ModelCIF templates and UniProt DBRefs, 1-byte AFDB audit revision presets, and sparse `-0.000` coordinate signs for exact mmCIF roundtrip fidelity.

---

## 2. Fixed Global Header (`BpdbV2Header`, 64 Bytes)

All multi-byte integer fields are little-endian. `sizeof(BpdbV2Header)` is statically asserted to be 64 bytes.

| Offset | Size | Type | Field | Description |
| :--- | ---: | :--- | :--- | :--- |
| `+0` | 4 | `char[4]` | `magic` | ASCII `"bPDB"` (`0x62 0x50 0x44 0x42`). |
| `+4` | 1 | `uint8_t` | `version` | Format version: `2`. |
| `+5` | 1 | `uint8_t` | `chain_id` | Single ASCII chain identifier (e.g., `'A'`). |
| `+6` | 1 | `uint8_t` | `c_term_oxt` | `1` if the C-terminal residue includes `OXT`; `0` otherwise. |
| `+7` | 1 | `uint8_t` | `cold_codec` | Cold Metadata compression codec: `0` = raw (uncompressed). |
| `+8` | 8 | `uint64_t` | `flags` | 64-bit format mode and cold-section presence bitmask (`FLAG_*`). |
| `+16` | 4 | `uint32_t` | `total_res` | Number of polymer residues ($N_{\text{res}}$). |
| `+20` | 4 | `uint32_t` | `total_atoms` | Total heavy atoms in the coordinate stream ($N_{\text{atoms}}$). |
| `+24` | 4 | `int32_t` | `first_seq_id` | Residue sequence ID of residue index 0 (typically `1`). |
| `+28` | 4 | `uint32_t` | `escape_count` | Number of 12-byte `EscapeAtom` entries ($N_{\text{esc}}$). |
| `+32` | 12 | `int32_t[3]` | `origin[3]` | Exact fixed-point coordinates of Atom 0 (`(x, y, z) * 1000`, milli-Å). |
| `+44` | 4 | `uint32_t` | `cold_offset` | Byte offset from start of stream where Cold Metadata begins. |
| `+48` | 12 | `char[12]` | `accession` | Null-padded UniProt accession (stripped of `AF-` prefix and `-F<n>` suffix). |
| `+60` | 2 | `uint16_t` | `frag_num` | AlphaFoldDB fragment index (e.g., `1` for `-F1`). |
| `+62` | 2 | `uint16_t` | `reserved` | Reserved; must be `0`. |

### 2.1 Global Flags Bitmask (`hdr.flags`)

| Bit | Mask | Constant | Meaning |
| ---: | :--- | :--- | :--- |
| 0 | `1ULL << 0` | `FLAG_AF_MODE` | Predicted structure (occupancy implicitly `1.00`, formal charge `?`). |
| 1 | `1ULL << 1` | `FLAG_AFDB` | Official AlphaFoldDB entry; enables compact AFDB cold sections and boilerplate mmCIF reconstruction. |
| 2 | `1ULL << 2` | `FLAG_HYDROGENS` | Explicit hydrogen atoms present (reserved for all-atom streams). |
| 3 | `1ULL << 3` | `FLAG_CIF` | Source/target text representation is mmCIF (ModelCIF). |
| 4 | `1ULL << 4` | `FLAG_AF2` | AlphaFold v2 semantics (e.g., `-0.000` sign preservation). |
| 5 | `1ULL << 5` | `FLAG_FRAGMENTED` | Entry belongs to a multi-fragment protein (`frag_num > 1`, or set by shard ingestor on `-F1` when sibling fragments exist). |
| 20 | `1ULL << 20` | `FLAG_RESOLUTION` | Cold tier contains 2-byte `uint16_t global_plddt` (scaled `×100`). |
| 34 | `1ULL << 34` | `FLAG_ENTITY` | Cold tier contains `_entity.pdbx_description`. |
| 48 | `1ULL << 48` | `FLAG_STRUCT` | Cold tier contains `_struct_conf` secondary structure records. |
| 49 | `1ULL << 49` | `FLAG_TEMPLATES` | Cold tier contains ModelCIF template TLV block. |
| 50 | `1ULL << 50` | `FLAG_DBREF` | Cold tier contains reference database (UniProt) TLV block. |
| 51 | `1ULL << 51` | `FLAG_AUDIT` | Cold tier contains `_pdbx_audit_revision_*` preset bitmask or raw history block. |
| 52 | `1ULL << 52` | `FLAG_NEG_ZERO` | Cold tier contains sparse `-0.000` coordinate sign table. |

---

## 3. Hot Geometry Tier (`[64 .. cold_offset)`)

Define 4-byte alignment as:

$$\operatorname{align4}(n) = (n + 3) \mathbin{\&} \sim 3$$

### 3.1 Sequence Array
* **Offset**: `64`
* **Size**: `align4(total_res)` bytes (trailing alignment bytes zero-filled).
* **Format**: 1 byte per residue holding the uppercase ASCII 1-letter amino acid code (`'A'`–`'Z'`). During AFDB parsing, `'U'` (selenocysteine) is mapped to `'C'` and `'O'` (pyrrolysine) is mapped to `'X'`.

### 3.2 Per-Residue pLDDT Array
* **Offset**: `64 + align4(total_res)`
* **Size**: `align4(total_res * 2)` bytes (trailing alignment bytes zero-filled).
* **Format**: Contiguous `uint16_t[total_res]` array of per-residue pLDDT scores scaled by `100` (e.g., `80.38` → `8038`), clamped to `[0, 65535]`. Every heavy atom in residue `r` shares `plddt[r]`.

### 3.3 12-Bit Bonded-Parent Delta Coordinate Stream
* **Offset**: `64 + align4(total_res) + align4(total_res * 2)`
* **Size**:
  $$\text{coord\_len} = \begin{cases} 0 & \text{if } N_{\text{atoms}} \le 1 \\ \operatorname{align4}\!\left(\left\lceil \frac{N_{\text{atoms}} - 1}{2} \right\rceil \cdot 9 + 7\right) & \text{if } N_{\text{atoms}} > 1 \end{cases}$$
  The trailing `+7` zero-padded bytes guarantee at least **16 readable bytes** from the start of the final 9-byte pair, allowing unmasked 128-bit SIMD loads (`_mm_loadu_si128` / `vld1q_u8`) without per-iteration bounds checks or buffer overreads.

#### Fixed-Point Representation & Topological Atom Order
* All Cartesian coordinates are quantized to exact **milli-Ångströms**: $q = \operatorname{round}(c \times 1000) \in \mathbb{Z}_{32}$.
* **Atom 0** (`N` of residue 0) is stored verbatim in `hdr.origin[3]`.
* The remaining $N_\Delta = N_{\text{atoms}} - 1$ atoms are emitted residue-by-residue in the **topological bonded-parent order** defined in `DICT[]` (`src/chem/dict.c`):
  * Backbone atoms always occupy indices `[0: N, 1: CA, 2: C, 3: O]`.
  * Side-chain heavy atoms follow in bonded-tree order such that every atom's covalently bonded parent has a strictly smaller index (`parent[a] < a`).
  * If `hdr.c_term_oxt == 1`, `OXT` is emitted as the final atom of residue `total_res - 1` with `parent = 2` (`C`).
  * For every residue $r \ge 1$, backbone `N` (index `0`) uses sentinel `BPDB_PARENT_PREV_C = 29`, referencing residue $r - 1$'s backbone `C` (index `2`).
  * When reconstructing mmCIF `_atom_site` records, atoms within each residue are re-emitted in standard `ATOM37_ORDER`.

#### Bonded-Parent Deltas & Escape Sentinel
For each atom $i \ge 1$, the encoder computes the signed milli-Ångström delta relative to its bonded parent:

$$\Delta = (x_i - x_{\text{parent}(i)},\; y_i - y_{\text{parent}(i)},\; z_i - z_{\text{parent}(i)})$$

* **Normal Delta**: $\Delta_x, \Delta_y, \Delta_z \in [-2047, +2047]$ ($\pm 2.047\text{ Å}$ per axis). Because covalent heavy-atom bonds span $\sim 1.2\text{–}1.85\text{ Å}$, $>99.95\%$ of AFDB structures require zero escapes.
* **Escape Sentinel**: If any axis lies outside $[-2047, +2047]$, the delta stream stores `(-2048, 0, 0)` (`0x800, 0x000, 0x000`) and appends the atom's exact absolute coordinates `{x, y, z}` to the **Escape Atom Table**.

#### 9-Byte Pair Bit-Packing (4.5 Bytes / Atom)
Consecutive deltas $\mathbf{d}_0 = (d_{0x}, d_{0y}, d_{0z})$ and $\mathbf{d}_1 = (d_{1x}, d_{1y}, d_{1z})$ — six 12-bit two's-complement signed integers (72 bits total) — are packed into 9 bytes:

```text
a0 (36 bits) = (d0.x & 0xFFF) | ((d0.y & 0xFFF) << 12) | ((d0.z & 0xFFF) << 24)
a1 (36 bits) = (d1.x & 0xFFF) | ((d1.y & 0xFFF) << 12) | ((d1.z & 0xFFF) << 24)

Byte:   [   0   ] [   1   ] [   2   ] [   3   ] [       4       ] [   5   ] [   6   ] [   7   ] [   8   ]
Bits:   a0[7:0]   a0[15:8]  a0 a0 a1[3:0]|a0 a1[11:4]  a1 a1 a1
```

If $N_\Delta$ is odd, $\mathbf{d}_1$ in the final pair is zero-filled `(0, 0, 0)`.

#### SIMD Unpacking Layout (`bpdb_decode_coords_f32`)
Each 9-byte pair is unpacked into two `int32x4` vectors (`[dx, dy, dz, 0]`) using a single 128-bit load at `p`:
1. **Byte shuffle** into 16-bit lanes using indices `{0,1, 1,2, 3,4, -1,-1,  4,5, 6,7, 7,8, -1,-1}`.
2. **Nibble alignment** (left-shift odd-nibble lanes by 4):
   * `x86-64 SSE4.1`: `_mm_mullo_epi16(shuffled, {16, 1, 16, 0,  1, 16, 1, 0})`
   * `AArch64 NEON`: `vshlq_s16(shuffled, {4, 0, 4, 0,  0, 4, 0, 0})`
3. **Sign extension**: Arithmetic right shift by 4 across all eight 16-bit lanes (`_mm_srai_epi16` / `vshrq_n_s16`).
4. **Widening**: Low 4 lanes → `int32x4` $\mathbf{d}_0$, high 4 lanes → `int32x4` $\mathbf{d}_1$.
5. **Fast Path**: When `hdr.escape_count == 0`, `bpdb_decode_coords_f32()` executes a branchless loop with unrolled backbone accumulation (`prev_C -> N -> CA -> C -> O`) and converts to `float4` via `vec * 0.001f`.

### 3.4 Escape Atom Table
* **Offset**: `coord_off + coord_len`
* **Size**: `escape_count * 12` bytes.
* **Format**: Contiguous array of `EscapeAtom { int32_t x, y, z; }` consumed sequentially whenever a decoded delta has `d[0] == -2048`.
* **Validation Invariant**: `coord_read_v2_hot()` verifies `magic == "bPDB"`, `version == 2`, `64 <= cold_offset <= buf_len`, and `esc_off + escape_count * 12 == hdr.cold_offset`.

---

## 4. Cold Metadata Tier (`[cold_offset .. EOF)`)

Cold metadata is tightly packed with no inter-field alignment padding. Wire strings are raw ASCII without a trailing `NUL` byte.

Several cold sections use **`varlen`** unsigned length/integer encoding:
* Values `0..254`: single `uint8_t` byte (`b < 0xFF`).
* Values `255..65535`: sentinel byte `0xFF` followed by little-endian `uint16_t`.

### 4.1 Column Formatting Widths (Always Present — 8 Bytes)
`uint8_t fmt_widths[8]` stores the exact whitespace-padded column widths detected from the source `_atom_site` loop for byte-exact mmCIF reconstruction:

| Index | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `_atom_site` Field | `id` | `label_atom_id` | `label_seq_id` | `Cartn_x` | `Cartn_y` | `Cartn_z` | `B_iso_or_equiv` | `pdbx_sifts_xref_db_num` |

### 4.2 Flag-Gated Metadata Sections (Written in Exact Order)

#### 1. `FLAG_RESOLUTION` (Bit 20)
* `uint16_t global_plddt` (2 bytes) — `_ma_qa_metric_global.metric_value` stored as fixed-point `percent * 100` (e.g., `93.67` → `9367`), consistent with per-residue pLDDT.

#### 2. `FLAG_ENTITY` (Bit 34)
* `varlen dlen` followed by `char[dlen]` — `_entity.pdbx_description`.
* **Compact sentinel**: `dlen == 0x00` decodes directly to `"Uncharacterized protein"`. Entries with no description omit `FLAG_ENTITY`.

#### 3. `FLAG_STRUCT` (Bit 48)
Secondary structure annotations (`_struct_conf`):
* **AFDB Compact Mode (`FLAG_AFDB` set)**:
  * `varlen count`
  * Repeated `count` times (typically 3 bytes per segment):
    * `varlen delta_start` — `start_idx - prev_start_idx` (`prev_start_idx` starts at `0`, where `start_idx = start.resSeq - first_seq_id`)
    * `varlen length` — `end.resSeq - start.resSeq + 1`
    * `uint8_t type_id` — encoded DSSP conformation ID
  * Start/end residue types and chain IDs are reconstructed directly from the Hot Sequence Array and `hdr.chain_id`.
* **Generic Mode (`FLAG_AFDB` clear)**:
  * `uint64_t count` followed by `MetaStructConf[count]` raw structs.
* **`type_id` Bit Layout**:
  * Category (`type_id & 0xE0`): `0x00` = Other/Helix (`0x00` = `OTHER`), `0x20` = Strand (`STRN`), `0x40` = Turn, `0x60` = Bend (`BEND`).
  * Subtype (`type_id & 0x1F`):
    * **Turn (`0x40`)**: `1` = `TURN_TY1_P`, `2` = `TURN_TY1P_P`, `3` = `TURN_TY2_P`, `4` = `TURN_TY2P_P`, `5` = `TURN_TY3_P`, `6` = `TURN_TY3P_P`, `7` = `TURN_P`.
    * **Helix (`0x00`)**: `(base << 1) | lh`, where `lh = 0` (right-handed `RH`) or `1` (left-handed `LH`), and `base` is `1`=`AL`, `2`=`OM`, `3`=`PI`, `4`=`GA`, `5`=`3T`, `9`=`27`, `10`=`PP`, `15`=`HELX_P`.

#### 4. `FLAG_TEMPLATES` (Bit 49) — Forward-Compatible TLV Block
* `uint32_t payload_bytes` — total byte length of everything following this field.
* `uint32_t template_count` — number of template records.
* Repeated `template_count` times: sequence of `(uint8_t tag, uint8_t len, uint8_t value[len])` entries terminated by `0x00` (`TPL_TAG_END`). Unknown tags are skipped by advancing `len` bytes. Common defaults are omitted on the wire:

| Tag | Constant | Value Encoding | Omission Default |
| :--- | :--- | :--- | :--- |
| `0x00` | `TPL_TAG_END` | 1-byte record terminator (no `len`) | — |
| `0x01` | `TPL_TAG_ID` | `uint16_t template_id` (`len = 2`) | `index + 1` (sequential) |
| `0x02` | `TPL_TAG_ACCESSION` | ASCII bytes (`len = 1..255`) | Empty / absent |
| `0x03` | `TPL_TAG_DB_NAME` | ASCII bytes (`len = 1..255`) | `"PDB"` |
| `0x04` | `TPL_TAG_ASYM_ID` | ASCII bytes (`len = 1..255`) | `"A"` |

#### 5. `FLAG_DBREF` (Bit 50) — Varlen TLV Block
* `uint32_t payload_bytes` — total byte length of everything following this field.
* Sequence of `(uint8_t tag, varlen len, uint8_t value[len])` entries terminated by `0x00` (`DBREF_TAG_END`). Unknown tags are skipped by advancing `len` bytes. Defaults/duplicates are omitted on the wire:

| Tag | Constant | Value Encoding | Omission Default |
| :--- | :--- | :--- | :--- |
| `0x00` | `DBREF_TAG_END` | 1-byte record terminator | — |
| `0x01` | `DBREF_TAG_DB_NAME` | ASCII string | `"UNP"` |
| `0x02` | `DBREF_TAG_ACCESSION` | ASCII string | `hdr.accession` |
| `0x03` | `DBREF_TAG_ID_CODE` | ASCII string (e.g., `"L2HDH_ACCRE"`) | Empty or encoded via `0x0E` |
| `0x04` | `DBREF_TAG_DB_ALIGN` | `int32_t[2]`: `db_seq_begin`, `db_seq_end` | `[first_seq_id, first_seq_id + total_res - 1]` |
| `0x05` | `DBREF_TAG_SEQ_ALIGN` | `int32_t[2]`: `seq_begin`, `seq_end` | `(0, 0)` |
| `0x06` | `DBREF_TAG_CHAIN` | `char` (1 byte) | `'A'` |
| `0x07` | `DBREF_TAG_ICODES` | `char[4]`: `iCode_beg/end`, `db_iCode_beg/end` | `"????"` |
| `0x08` | `DBREF_TAG_CHECKSUM` | `uint64_t` (64-bit UniProt sequence CRC) | `0` |
| `0x09` | `DBREF_TAG_TAX_ID` | `uint32_t` (NCBI taxonomy ID) | `0` |
| `0x0A` | `DBREF_TAG_VERSION_DATE` | 3-byte packed `{uint8_t yr-1973, mo, day}` (or raw ASCII if `len != 3`) | Empty / absent |
| `0x0B` | `DBREF_TAG_GENE` | ASCII string | Empty / absent |
| `0x0C` | `DBREF_TAG_ORGANISM` | ASCII string (stored verbatim; tax ID is not bijective) | Empty / absent |
| `0x0D` | `DBREF_TAG_ISOFORM` | ASCII string | Absent or `"?"` |
| `0x0E` | `DBREF_TAG_ID_SUFFIX` | ASCII suffix reconstructed as `<accession>_<suffix>` (TrEMBL ID compaction) | Used when `db_id_code` has `<accession>_` prefix |

#### 6. `FLAG_AUDIT` (Bit 51)
Starts with a single `uint8_t audit_mask`:
* **Preset Bitmask Mode (`0x01..0x3F`)**: 1 byte total; each bit selects a canonical AFDB release row from `AFDB_AUDIT_PRESETS[]` (`src/cif/afdb_strings.c`):
  * Bit `0` (`0x01` `AUDIT_MASK_V1`): `2021-07-01` v1.0 Initial release
  * Bit `1` (`0x02` `AUDIT_MASK_V2`): `2021-12-09` v2.0 Remediation (SwissProt)
  * Bit `2` (`0x04` `AUDIT_MASK_V3`): `2022-06-01` v3.0 Remediation (UniProt)
  * Bit `3` (`0x08` `AUDIT_MASK_V4`): `2022-09-30` v4.0 Remediation (accuracy)
  * Bit `4` (`0x10` `AUDIT_MASK_V5`): `2025-03-31` v5.0 Release (UniProt 2024_06)
  * Bit `5` (`0x20` `AUDIT_MASK_V6`): `2025-08-01` v6.0 Release (UniProt 2025_03)
* **Raw Inline Fallback (`0xFF` `AUDIT_MASK_RAW`)**: Followed by `uint32_t revision_count` and `revision_count` records of:
  * Packed 8-byte `MetaAuditRevision`: `{ uint8_t year (YYYY - 1973), month, day, major, minor, type_len; uint16_t desc_len; }`
  * `char[desc_len]` description, followed by `char[type_len]` revision type.

#### 7. `FLAG_NEG_ZERO` (Bit 52)
* `uint32_t neg_zero_count`
* Sorted `uint32_t[neg_zero_count]` entries packed as `(stream_atom_idx << 2) | axis`, where `stream_atom_idx` is the 0-based atom index in topological stream order and `axis` is `1` (`X`), `2` (`Y`), or `3` (`Z`). Preserves `-0.000` vs `0.000` during mmCIF roundtrip.

---

## 5. Public C API (`include/bpdb.h`)

```c
// Lossless mmCIF <-> bPDB v2
size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf);
void   decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f);

// Reusable thread-local encoder context (zero per-file scratch reallocation in batch pipelines)
typedef struct BpdbEncoderCtx BpdbEncoderCtx;
BpdbEncoderCtx* bpdb_encoder_create(void);
void            bpdb_encoder_destroy(BpdbEncoderCtx* ctx);
size_t          bpdb_encode_reuse(BpdbEncoderCtx* ctx, const char* cif_data, size_t cif_len,
                                  const uint8_t** out_ptr);

// Zero-allocation Hot Geometry API (mmap-friendly)
bool bpdb_open_hot(const uint8_t* buf, size_t buf_len, BpdbHotView* out_view);
bool bpdb_decode_coords_f32(const uint8_t* buf, size_t buf_len, float* out_xyzw);
```

---

## 6. Building & Testing

```bash
# Build libbpdb.a and run lossless roundtrip test against the fixture
make test

# Build the streaming AFDB tar.gz batch harness (recommended: libdeflate)
make harness USE_LIBDEFLATE=1

# Clean build artifacts and temporary outputs
make clean
```

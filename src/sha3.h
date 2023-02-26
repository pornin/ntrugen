#ifndef SHA3_H__
#define SHA3_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Context for a SHAKE computation. Contents are opaque.
 * Contents are pure data with no pointer; they need not be released
 * explicitly and don't reference any other allocated resource. The
 * caller is responsible for allocating the context structure itself,
 * typically on the stack.
 * A running state can be cloned by copying the structure; this is
 * useful if "partial hashes" (hash of data processed so far) are
 * needed, without preventing injecting extra bytes later on.
 */
typedef struct {
	uint64_t A[25];
	unsigned dptr, rate;
} shake_context;

/*
 * Initialize a SHAKE context to its initial state. The state is
 * then ready to receive data (with shake_inject()).
 *
 * The "size" parameter should be 128 for SHAKE128, 256 for SHAKE256.
 * This is half of the internal parameter known as "capacity" (SHAKE128
 * works on an internal 256-bit capacity, SHAKE256 uses a 512-bit
 * capacity).
 */
void shake_init(shake_context *sc, unsigned size);

/*
 * Inject some data bytes into the SHAKE context ("absorb" operation).
 * This function can be called several times, to inject several chunks
 * of data of arbitrary length.
 */
void shake_inject(shake_context *sc, const void *data, size_t len);

/*
 * Flip the SHAKE state to output mode. After this call, shake_inject()
 * can no longer be called on the context, but shake_extract() can be
 * called.
 *
 * Flipping is one-way; a given context can be converted back to input
 * mode only by initializing it again, which forgets all previously
 * injected data.
 */
void shake_flip(shake_context *sc);

/*
 * Extract bytes from the SHAKE context ("squeeze" operation). The
 * context must have been flipped to output mode (with shake_flip()).
 * Arbitrary amounts of data can be extracted, in one or several calls
 * to this function.
 */
void shake_extract(shake_context *sc, void *out, size_t len);

/*
 * Structure holding four SHAKE context simultaneously. This is used for
 * the "x4" mode, which outputs 64-bit words upon extraction, and may
 * leverage architecture-specific vector instruction for better performance.
 */
typedef struct {
	uint64_t A[100];
	unsigned dptr, rate;
} shake_x4_context;

/*
 * Flip four input SHAKE context to output mode. The four provided contexts
 * are unmodified, but a new shake_x4_context structure is filled with
 * their "flipped" counterparts. That structure may be used with
 * shake_x4_extract_words() to obtain the four outputs as 64-bit words
 * (the outputs are interleaved with a 64-bit granularity).
 *
 * All four input SHAKE context MUST use the same internal capacity
 * (i.e. do not mix SHAKE128 and SHAKE256 contexts). The internal
 * capacity MUST be a multiple of 64 bits.
 */
void shake_x4_flip(shake_x4_context *scx4, const shake_context *sc_in);

/*
 * Obtain num_x4 groups of four 64-bit words from the provided parallel
 * SHAKE context. Each group contains one word from each of the four
 * internal SHAKE instances, in the order they were provided when
 * calling shake_x4_flip(). Each 64-bit word corresponds to the little
 * endian interpretation of corresponding 8-byte output.
 */
void shake_x4_extract_words(shake_x4_context *scx4,
	uint64_t *dst, size_t num_x4);

/*
 * Context for SHA3 computations. Contents are opaque.
 * A running state can be cloned by copying the structure; this is
 * useful if "partial hashes" (hash of data processed so far) are
 * needed, without preventing injecting extra bytes later on.
 */
typedef shake_context sha3_context;

/*
 * Initialize a SHA3 context, for a given output size (in bits), e.g.
 * set size to 256 for SHA3-256.
 */
void sha3_init(sha3_context *sc, unsigned size);

/*
 * Update a SHA3 context with some bytes.
 */
void sha3_update(sha3_context *sc, const void *in, size_t len);

/*
 * Finalize a SHA3 computation. The hash output is written in dst[],
 * with a size that depends on the one provided when the context was
 * last initialized.
 *
 * The context is modified. If a new hash must be computed, the context
 * must first be reinitialized explicitly.
 */
void sha3_close(sha3_context *sc, void *out);

#ifdef __cplusplus
}
#endif

#endif

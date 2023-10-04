# BRRR

`cargo r -r -- file.fasta`

## Developper's notes

Minimizers are computed using a monotone queue, the order is based on a hash function which can be seeded using the `new_with_seed` method.

The `bloom` module provides an implementation of Bloom filters and cascading Bloom filters.
These Bloom filters compute the hashes based on two hash functions using double hashing, once again the hash functions can be seeded using the `new_with_seed` method.
In order to improve cache-efficiency, the hashes associated to an element are all mapped to a single block that fits in cache.

The `dashbloom` module provides a drop-in replacement of the (cascading) Bloom filters from the `bloom` module, and allows concurrent access to the filters by different threads.
This is done by using the first few bits of the hashes to dispatch the elements between smaller Bloom filters that are thread-safe.

In order to process the reads with multiple threads, just replace the `reads.process()` function by `reads.parallel_process()` while specifying the number of threads and the size of the queue.

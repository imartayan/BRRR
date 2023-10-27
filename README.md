# Brutal Read Rewrite in Rust

## Usage

```md
cargo r -r -- [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file (.fasta, .fa)

Options:
  -o, --output <OUTPUT>        Output file (defaults to <input>.cor.<ext>)
  -t, --threads <THREADS>      Number of threads (defaults to all available threads)
  -m, --memory <MEMORY>        Memory (in MB) allocated to Bloom filters (defaults to input size)
  -a, --abundance <ABUNDANCE>  Abundance above which k-mers are solid [default: 5]
  -H, --hashes <HASHES>        Number of hashes used in Bloom filters [default: 3]
  -s, --seed <SEED>            Seed used for hash functions [default: 101010]
  -h, --help                   Print help
  -V, --version                Print version
```

By default `K=31` and `M=21` are fixed, but you can specify other values as follows:
```sh
K=15 M=7 cargo r -r -- [OPTIONS] <INPUT>
```

## Developer's notes

Minimizers are computed using a monotone queue (with lookup in *O(1)* and insertion in amortized *O(1)*), the order is based on a hash function which can be seeded using `-s`.

The `bloom` module provides an implementation of Bloom filters, cascading Bloom filters and counting Bloom filters.
These Bloom filters compute the hashes based on two hash functions (also seeded with `-s`) using double hashing.
In order to improve cache-efficiency, the hashes associated to an element are all mapped to a single block that fits in cache.

The `dashbloom` module provides a drop-in replacement of the different kinds of Bloom filters implemented in `bloom`, and allows concurrent access to the filters by different threads.
This is done by using the first few bits of the hashes to dispatch the elements between smaller Bloom filters that are thread-safe.

In order to process the reads with multiple threads, simply replace the `reads.process` function by `reads.parallel_process` while specifying the number of threads and the size of the queue.
If the threads need to send a result (such as a buffer when rewriting reads), `parallel_process_result` should be used instead.

Note that rewriting the reads using multiple threads may not preserve the original order of the reads.

The main strategy used in the correction method uses a bidirectional BFS (implemented in `find_path`) to find a solid path between the two closest solid k-mers.

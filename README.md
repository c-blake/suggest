This is a from-scratch implementation in Nim of Wolf Garbe's Symmetric Delete
Algorithm for correct spelling suggestions.  The basic idea is simple: corpus
words within edit distance N of a query can only be *at most* N units longer or
N units shorter and the shorter must be derived from deletes of longer strings.
So, build a map from all shortenings to all corpus words which generate them.
When queried with a string, we can then lookup all possible edits lengthening
the query into a corpus word.  We also on-the-fly compute all shortenings of a
query (& shortenings of lengthenings).  From both sets, we filter "maybe within
N edits of a corpus word" to "actually within N edits".  The filter can be any
"distance" successfully bounded by N-indels.  This idea is like Monte Carlo
numerical integration of shapes within easy bounding boxes (but this is
deterministic & points which pass are reported, not just counted).

While this does allow for fast queries, it takes takes a long time and a lot of
space to make a lengthenings table compared to *thousands* of linear scans of a
large-ish correct word corpus.  So, it is a useful strategy if you A) can save
the big map to disk **and very** efficiently load it **and/or** B) have **many**
queries to amortize build costs over.  Rebuilding is lame while "real" DB query
latency is hostile.  So, this module does an efficient external file format with
five files to "mmap & go":
```
A .tabl pointing to (keyIx->varlen.keys, .sugg=varlen[array[CNo]])
and a .meta(ix,cnt) file pointing to varlen .corp.
```
`varlen[arr[CNo]]` is an typical allocation arena with early entries the heads
of per-list-size free lists.

I originally wrote this to understand and perhaps debunk SymSpell, though the
work has (sort of) validated it.  Or not.  Let you, dear reader, are the judge.
I think that I have at least found information I didn't see elsewhere that
brackets its applicability which is worth letting people know about.  In
particular, SymSpell offers only modest speed-up vs-linear scan at large
(4,5,..) edit distances of a medium- sized (40 kWord) corpus, as shown in ![this
plot.](https://raw.githubusercontent.com/c-blake/suggest/master/scanVsymspellD5.png)
This does roughly contradict Garbe's "large distance, large dictionary" sales
pitch.  False positive rates for d>3 probably makes that regime uninteresting.
Still, SymSpell benefit remains only 3.5x-ish for 80 kWord which is not great.
Indeed, multi-core storage/processing optimizations on both linear scan and
symspell querying might even nullify such a small advantage.

The basic experimental set up is to use "frequency\_dictionary\_en\_82\_765.txt"
from the SymSpell repository as our input.  We create synthetic "batches of
typos" by sampling words from that very same distribution and then deleting one
or more random characters.  Create 2000 synthetic batches and feed these into a
'compare' mode of 'suggest' that times both a linear scan approach and an mmap &
query SymSpell approach.  The averages should thus be pretty representative.
Batch size of 6 typos per document/user interaction felt about right, but typo
rates probably vary a lot.  The script is provided for users wanting to study
how things vary.  The error bars are the standard deviation of the mean.  The
distribution is wide with 95th percentile times often 4X the 5th percentile.

Both linear scanning and symspell querying have an additional optimization of
shrinking the max distance passed to optimized edit distance calculators once
"enough" correct suggestions have been found at lesser distances.  This can
speed up such distance computations somewhat, especially for linear scans of
short words near popular portions of word-space or when very few matches are
requested.

You can see the general scaling of SymSpell costs with max distance and
vocabulary size from ![this
plot.](https://raw.githubusercontent.com/c-blake/suggest/master/scanVsymspell4k.png)

Along the way, I also found that SymSpell is a pretty good stress test for a
some system-layer functionality - in particular memory allocators, string hash
functions, and even virtual memory operation.  Indeed, unless you really need
that last 10x-50x speed boost from 100s of microseconds to 10s of microseconds,
it is likely not advisable to try to implement SymSpell (at least without the
systems layer warnings in this README).  A simple linear scan of a corpus which
has been pre-sorted by descending recommendation order is fast enough and very
simple - you needn't even sort results.  You simply bin them by edit distance.

For the curious about the systems layer problems, sensitivity to string hash
functions is perhaps the most obvious.  Because each key in the main table is
formed by 1-deletions (and more) of other keys the amount of key similarity is
an almost hash-attack level end-point.  One doesn't quite need cryptographically
secure hash functions, but almost any "weak but fast" will be defeated.  The
symptom is the usual one - long collision chains and slow operations on what is
usually a gigantic table.  My timing results use a variant of the late Robert
Uzgalis' BuzHash that runs in (1.75 cycles + 0.63 cycles/byte).  The default Nim
hash also works all right, but will have not quite as good numbers.  (For one,
it's 4.3cyc+4.2 cyc/byte, but also not quite as randomizing.)  I tried hard to
have SymSpell be cast in the best light possible here, and that meant using the
best hash.  Similarly, I did not try to implement the prefix optimization that
https://github.com/wolfgarbe/symspell discusses which slows things down a little
at some significant memory savings.

The memory allocator problem comes from a long-tailed distribution of how many
"correct suggestions" any given typo has.  The distribution shape defeats most
attempts to "speed-up" memory allocators with, say, "power of two" spaced region
sizes.  Almost any spacing besides the minimal one results in very low space
utilization by suggestion lists.  Thankfully, minimal spacing is fast enough.
Indeed, an early non-persistent version of the code blew up most GCs Nim offers.
Only the Boehm-Demers-Weiser garbage collector actually allowing completion of a
table build in reasonable time.

Finally, while 4096 byte virtual memory pages are rarely a performance obstacle,
the size and access pattern of symspell queries is particularly hostile to use
from a fresh mmap.  For larger dictionaries and distances <=~ 3, using 2M pages
so-called "huge TLB" pages resulted in >2x speed-ups for a "fresh mapping", as
can be seen in ![this
graph](https://raw.githubusercontent.com/c-blake/suggest/master/4kVs2M.png)
That relative speed-up of large pages does owe to the small absolute time
SymSpell queries take, of course, but is still indicative of thousands of
non-local 4k page accesses which will become relevant in later discussions.

It also bears mentioning that table building time is still costly in this fairly
optimized implementation.  When using a Linux tmpfs RAM filesystem the resource
usage looks like this:
```
maxLevD dt(sec) size(bytes) Misspells
      1   0.182  19,456,529    647786
      2   0.957  76,853,412   2562636
      3   3.203 175,542,825   6852017
      4   8.127 364,447,609  13981420
      5  17.034 694,357,352  23487763
```
Those times could be sped up 1.2-1.3x by clairvoyantly pre-sizing the table to
avoid hash table resizes.  Some reasonably precise formula for the number of
delete misspellings for a given (corpus, distance) might be able to provide this
modest speed-up systematically.  It is also very unlikely that "less persistent"
implementations can achieve build times much better.
```
maxLevD dt(sec) size(bytes) Misspells
      1   0.141  19,456,529    647786
      2   0.710  76,853,412   2562636
      3   2.589 175,542,825   6852017
      4   6.654 364,447,609  13981420
      5  13.891 694,357,352  23487763
```
Also Garbe claims significant space reduction for modest slow downs which should
also help build times.  The main take-away though, is just making concrete my
"long time and a lot of space" from the introductory paragraphs.  This puts real
pressure on saving the answer of this build, especially for larger dictionaries
with longer words.  One really does need thousands of future queries before the
investment in build time pays off in query performance.

In comparison, the corpus file alone is a mere 752,702 bytes and can be scanned
in hundreds of microseconds off a modern NVMe storage even *fully cold-cache*.
Cold-cache, non-RAM FS times for SymSpell, meanwhile degrade dramatically,
especially on high latency storage like Winchester drives.  Every corrections
query involves at least one, but sometimes hundreds or even thousands of random
accesses to the `.tabl` file for weakly related keys.  Arranging for locality to
such accesses seems quite challenging if not impossible.

For example, Cold-cache scanning on a 10 ms + 1e-5 ms/MB (aka 10 ms, 100 MB/s)
Winchester drive or network storage unit would probably be approximately
10 + 1e-5\*750e3 or about 17.5 ms, generally half the time of a single frame
of a typical video file.  Getting suggestions for multiple words can easily
share that IO as well, making amortized per-batch-of-6-typos more like 3 ms per
word.  Meanwhile, just one cold-cache suggestion on freshly opened data files
could take 1000s of random accesses or 10s of seconds (at 10 ms latencies).

TODO: I should instrument the code to print out random accesses per query again
and put some graph of those crazy distributions for d=1-5 here.

Multiple words also do not share the random aspects of SymSpell IO well and a
minute or multiple minutes is conceivable.  Indeed, if multiple words are
expected, the wisest IO strategy would be to get the whole data set paged into
RAM in a streaming sweep and then run the queries which would only take 7 sec
for the largest d=5 case.  This cold-cache scenario is yet another "system layer
performance fragility" of SymSpell, but notice that linear scan's 10s of ms can
now be 1000x faster than SymSpell's 10s of seconds.  Faster than 10 ms + 1e-5
ms/MB storage is becoming more common these days, but even so..to keep SymSpell
a performance winner in deployment one wants to ensure cached pages.

Of course, the above commentary also applies to *non-persistent* SymSpell
implementations in execution environments where swap/page files are possible.
Competition for memory may be beyond a developers control, and pre-paging
non-persistent data can be much harder than "cat spell.\* > /dev/null".
Also, SymSpell's large memory requirements make it probably one of the more
fierce memory competitors.

The TL;DR?  While a well-implemented SymSpell can indeed be faster than a
similarly well-implemented linear scan, it is far more "performance risky".
It could be 100x faster than a linear scan in some hot cache circumstances
or 1000x slower in other cold-cache circumstances.  Meanwhile, a cold-cache
linear scan might be only about 20x worse than hot-cache.  In these terms,
SymSpell is 5000x more performance risky than a linear scan.  All these risks
can be addressed, but the developer needs to be mindful of them.

The TL;DR;DR?  "YMMV from hell".  ;-)

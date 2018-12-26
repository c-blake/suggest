## A from-scratch re-implementation of Wolf Garbe's Symmetric Delete Algorithm
## for correct spelling suggestions.  The basic idea is simple: corpus words
## within edit distance N of a query can only be *at most* N units longer or N
## units shorter and the shorter must be derived from deletes of longer strings.
## So, build a map from all shortenings to all corpus words which generate them.
## When queried with a string, we can then lookup all possible edits lengthening
## the query into a corpus word.  We also on-the-fly compute all shortenings of
## a query (& shortenings of lengthenings).  From both sets, we filter "maybe
## within N edits of a corpus word" to "actually within N edits".  The filter
## can be any "distance" successfully bounded by N-indels.  This idea is like
## Monte Carlo numerical integration of shapes within easy bounding boxes (but
## this is deterministic & points which pass are reported, not just counted).
##
## While this allows for fast queries, it takes takes a long time and much space
## to make a lengthenings table.  So, it is a useful strategy if you A) can save
## the big map to disk AND VERY efficiently load it AND/OR B) have MANY queries
## to amortize build costs over.  Rebuilding is lame & "real" DB query latency
## is hostile.  So, this module does an efficient external file format w/5 files
## to mmap&go: A .tabl pointing to (keyIx->varlen.keys,.sugg=varlen[array[CNo]])
## and a .meta(ix,cnt) file pointing to varlen .corp.  varlen[arr[CNo]] is an
## allocation arena with early entries the heads of per-listSz free lists.

import hashes,tables,sets, os, times, memfiles, strutils, algorithm, math,random
type
  MyersPattern*[W] = object       ## Output of myersCompile(var Pattern[W],..)
    m: W                          ## Length of pattern string
    pm: seq[W]                    ## Match vector: pm[ch]&(1<<i)<->A[i]==ch

proc myersCompile*[T,W](pat: openArray[T], dummy: W): MyersPattern[W] =
  ## Compile pattern for either levenshtein or optimStrAlign.
  result.m  = pat.len
  result.pm = newSeq[W](int(T.high))
  for i in 0 ..< pat.len:
    result.pm[int(pat[i])] = result.pm[int(pat[i])] or (1 shl i)

proc levenshtein*[T,W](p: MyersPattern[W], t: openArray[T], maxDist: W): W =
  ## Myers BitVec Algo after Hyyro 2003 Paper version (hP..or 1) with an early
  ## len check against maxDist.
  if p.m == 0: return t.len             #Degenerate case
  let m = p.m
  let n = t.len
  if n-m >= maxDist or m-n >= maxDist:  #Paper tests > & returns +1, but..
    return maxDist                      #..want convention of other procs.
  let mxPBit = 1.W shl (m - 1)          #Bit representing change in last row
  var vP = -1.W                         #Vert change is initially all pos
  var vN = 0.W                          #I.e., no neg.
  var d  = m.W
  for j in 0 ..< n :
    let pmj = p.pm[int(t[j])]
    var d0 = (((pmj and vP) + vP) xor vP) or pmj or vN
    var hP = vN or not (d0 or vP)
    let hN = d0 and vP
    if   (hP and mxPBit) != 0: inc(d)
    elif (hN and mxPBit) != 0: dec(d)
    hP = (hP shl 1) or 1
    vP = (hN shl 1) or not(d0 or hP)
    vN = d0 and hP
  return min(maxDist, d)

proc optimStrAlign*[T,W](p: MyersPattern[W], t: openArray[T], maxDist: W): W =
  ## Myers BitVec Algo after Hyyro 2003 Paper version with early len ck adapted
  ## to optimStrAlign.
  if p.m == 0: return t.len             #Degenerate case
  let m = p.m
  let n = t.len
  if n-m >= maxDist or m-n >= maxDist:  #Paper tests > & returns +1, but..
    return maxDist                      #..want convention of other procs.
  let mxPBit = 1.W shl (m - 1)          #Bit representing change in last row
  var vP   = -1.W                       #Vert change is initially all pos
  var vN   = 0.W                        #I.e., no neg.
  var d    = m.W
  var d0   = vN
  var pmj0 = vN
  for j in 0 ..< n :
    let pmj = p.pm[int(t[j])]
    d0 = ((((not d0) and pmj) shl 1) and pmj0) or
           (((pmj and vP) + vP) xor vP) or pmj or vN
    pmj0 = pmj
    var hP = vN or not (d0 or vP)
    let hN = d0 and vP
    if   (hP and mxPBit) != 0: inc(d)
    elif (hN and mxPBit) != 0: dec(d)
    hP = (hP shl 1) or 1
    vP = (hN shl 1) or not(d0 or hP)
    vN = d0 and hP
  return min(maxDist, d)

const WHI* = 31                   ## Last byte = chars actually used by word
type
  Ix* = uint32                    ## integer type for byte offsets into files
  CNo* = uint32                   ## integer type for record offsets into corp
  Count* = float64                ## float32 ranks word prob fine; breaks@>2**24
  Word* = object
    n*: uint8
    d*: array[WHI, char]          ## Bounded array of *up to* WHI chars
  TabEnt* = object {.packed.}     ## Size=13 (should be corp:1 => 12; Nim bug?)
    keyI*: Ix                     ## Possible typo -> suggestion list
    keyN* {.bitsize:  8.}: uint8  ## typo len; (Could be hash w/len in .keys)
    corp* {.bitsize:  1.}: uint8  ## 1-bit flag indicating a real corpus word
    szIx* {.bitsize:  7.}: uint8  ## Index into `sizes` (alloc);
    len*  {.bitsize: 16.}: uint16 ## Length of sugg list (used)
    sugg*: Ix                     ## pointer into file of suggestions lists
  SuggTab* = object {.packed.}    ## 32B header & table
    maxDist*: uint8               ## Within `maxDist` of corpus words
    unique* : CNo                 ## Track unique words add()d
    total*  : float64             ## Track all words add()d
    len*    : Ix                  ## Non-empty entries in table
    tab* : UncheckedArray[TabEnt] ## Hash structured; 12B per entry
  Suggs* = UncheckedArray[CNo]    ## Suggs Arena; 1st 16 are heads of free lists
  Meta* = object {.packed.}
    cix*: Ix                      ## Byte index into .corp for corpus word
    cnt*: Count                   ## Fixed size statistics on the corpus word
  Metas* = UncheckedArray[Meta]   ## Meta data for corpus words
  Suggestor* = object
    mode* : FileMode              ## Mode suggestor data files are in
    tabf* : MemFile               ## Backing store handle for table
    tabSz*: int                   ## Table size in TabEnt units; Calc from .size
    table*: ptr SuggTab           ## Run-time table handle
    keyf* : MemFile               ## TabEnt.keyI points to store for table keys.
    zKeys*: Ix                    ## USED bytes in keys file; alloc=keyf.size
    sugf* : MemFile               ## Backing store handle for words
    sugSz*: Ix                    ## suggs size in CNo units; Calc from .size
    suggs*: ptr Suggs             ## Run-time sugg list file handle
    metf* : MemFile               ## Backing store for (byteOffset,Count)
    metas*: ptr Metas             ## Run-time medata file handle
    nMeta*: CNo                   ## USED recs; alloc=metf.size div sizeof(Meta)
    corf* : MemFile               ## Backing store for (Count,CorpWord) @CNo's.
    zCorp*: Ix                    ## USED BYTES in corp file; alloc=corf.size
    saved*: array[5, int]         ## HugeTLBfs on Linux rounds st_size.  So, we
                                  ##..get via refr, but munmap wants .saved.
  Results* = seq[seq[CNo]]        ## Per-distance-from-query seq[corpusIx]
  DistanceKind* = enum lev, osa   ## Levenshtein, Optimal String Alignment
  ucArrCh* = ptr UncheckedArray[char]
const szTabEnt = 12 #sizeof(TabEnt) #XXX =13 instead of same as C sizeof==12.

var totDists = 0
proc distance*[W](p: MyersPattern[W], t: ptr char, n: int,
                  maxDist: W, kind=osa): W {.inline.} =
  inc(totDists)
  if kind == lev: p.levenshtein(toOpenArray(cast[ucArrCh](t), 0, n-1), maxDist)
  else        : p.optimStrAlign(toOpenArray(cast[ucArrCh](t), 0, n-1), maxDist)

proc `+%`(p: pointer, i: uint): pointer {.inline.} =  #Pointer arithmetic
  cast[pointer](cast[uint](p) + cast[uint](i))

template te[T](i: T): auto = s.table[].tab[int(i)]    #TableEnt accessor

#NOTE: .keyI("") != 0 because very first key added is always a real corpus word.
proc empty(e: TabEnt): bool {.inline.} = e.keyI == 0 and e.keyN == 0

proc toWord*(s: string): Word {.inline.} =
  result.n = uint8(min(WHI, s.len))
  copyMem(addr result.d[0], unsafeAddr s[0], result.n)

proc `$`*(w: Word): string {.inline.} =
  let n = int(w.n)
  result.setLen(n)
  if n > 0: copyMem(addr result[0], unsafeAddr w.d[0], n)

#PERSISTENT MMAP()d HASH TABLE
proc keyData*(s: Suggestor; i: Ix): pointer {.inline.} =
  if te(i).corp == 0:                   #Not corpus word => from .keys
    s.keyf.mem +% uint(te(i).keyI)
  else:                                 #corp=>metas->corp; +1 skips .n field
    s.corf.mem +% uint(s.metas[int(te(i).keyI)].cix + 1)

proc toWord*(s: Suggestor; i: Ix): Word {.inline.} =
  result.n = te(i).keyN
  copyMem(addr result.d[0], s.keyData(i), result.n)

proc mcmp(a,b: pointer; n: csize): cint {.importc:"memcmp",header:"<string.h>".}
var cMax = 0
var tabFinds = 0
proc find*(s: Suggestor, w: Word): int =
  inc(tabFinds)
  let mask = s.tabSz - 1                #Vanilla linear probe hash search/insert
  var i = hash(w.d) and mask            #Initial probe
  let n = int(w.n)                      #Length of query word
  var c = 0
  while not te(i).empty:
    inc(c)
    if int(te(i).keyN)==n and mcmp(s.keyData(i.Ix), unsafeAddr w.d[0], n)==0:
      return i                          #Len & *only* then bytes equal => Found
    i = (i + 1) and mask                #The linear part of linear probing
  if c >= 64*10 and c > cMax:   #64 giant lg Sz; 10 giant statistical factor.
    cMax = c                            #Print @new max to not do too many msgs
    stderr.write "weak hash function: ", c, " loops to find empty for ", w, "\n"
  return -i - 1                         #Not Found, return -(insertion point)-1

proc tBytes*(nSlot:int): int {.inline.} = sizeof(SuggTab) + (nSlot * szTabEnt)
proc tSlots*(nByte:int): int {.inline.} = (nByte-sizeof(SuggTab)) div szTabEnt

proc growTab*(s: var Suggestor, w: Word): int =
  let oldSz = s.tabSz   #Triple tab[] space, copy old to top, zero new. Might..
  s.tabSz *= 2          #..avoid copy&zero w/more syscalls for marginal speedup.
  s.tabf.resize(tBytes(s.tabSz + oldSz))
  s.table = cast[ptr SuggTab](s.tabf.mem)
  copyMem(addr te(s.tabSz), addr te(0), oldSz * szTabEnt)
  zeroMem(addr te(0), s.tabSz * szTabEnt)
  for i in s.tabSz ..< s.tabSz + oldSz:             #For each possible OLD slot
    if te(i).empty: continue                        #skip empty
    let J = -s.find(s.toWord(Ix(i))) - 1            #find spot & copy ent to new
    let j = if J < 0: -J - 1 else: J            #Can falsePos for "" in NEW tab
    copyMem(addr te(j), addr te(i), szTabEnt)   #but cpOver of old non-empty wks
  s.tabf.resize(tBytes(s.tabSz))                    #Now shrink file & memory
  return -s.find(w) - 1                   #Return insertion spot in grown table

proc addKey*(s: var Suggestor, mim1: int, w: Word, c: bool): int =
  var i = -mim1 - 1                                         #Add key w/no suggs
  if int(s.table[].len) + 1 > 7 * (s.tabSz shr 3):          #> 87.5% => double
    i = s.growTab(w)
  s.table[].len += 1
  if int(s.zKeys + Ix(w.n)) > s.keyf.size:          #Maybe grow .keys file
    s.keyf.resize(int((s.zKeys + Ix(w.n)) shl 1))
  copyMem(s.keyf.mem +% uint(s.zKeys), unsafeAddr w.d[0], w.n)
  te(i).keyI = s.zKeys                              #Point table at new key
  te(i).keyN = w.n
  s.zKeys += Ix(w.n)
  te(i).corp = if c: 1 else: 0
  return i

#ALLOCATION ARENA FOR PERSISTENT SUGGESTION LISTS
const  sizes = [ 0, 1,    2,    3,   4,   5,   6,   7,   8,   9,   10,  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 48, 56, 64, 85, 102, 113, 128, 146, 170, 204, 256, 341, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65535 ]
#const grows = [ 0, 1024, 1024, 682, 512, 409, 341, 292, 128, 113, 102, 92, 85, 78, 73, 68, 64, 60, 56, 53, 51, 48, 46, 44, 42, 40, 39, 37, 36, 35, 34, 33, 32, 25, 21, 18, 16, 12, 10,  9,   8,   7,   6,   5,   4,   3,   2,   1,    1,    1,    1,    1,     1,     1     ]
#Below is needed for 70%-ish utilization on the .sugg file. Could ditch grows[]
const  grows = [ 0, 1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,     1,     1     ]

proc growArea*(s: var Suggestor, head: uint16) =
  let n = Ix(grows[head] * sizes[head])         #How much to grow by; CNo units
  let oldSz = s.sugSz                           #oldSz; Ix(new data); CNo units
  s.sugSz += n                                  #Update min needed file length
  if int(s.sugSz) * sizeof(CNo) > s.sugf.size:  #Need more file space:
    s.sugf.resize((int(s.sugSz) * sizeof(CNo)) shl 1) #Double, remap and fix..
    s.suggs = cast[ptr Suggs](s.sugf.mem)             #..convenience field.
  s.suggs[head] = CNo(oldSz)                    #head = oldSz & thread rest
  if grows[head] >= 2:                          #New space from OS is all 0
    var ix = oldSz                              #NOTE: This free-list threading
    while ix < oldSz + Ix(grows[head] - 2):     #..is inactive with grows[]=1.
      s.suggs[ix] = CNo(ix + Ix(sizes[head]))
      ix += Ix(sizes[head])

proc free*(s: var Suggestor; head: uint16, i: Ix) {.inline.} =
  zeroMem(addr s.suggs[i], sizeof(Ix)*sizes[head])  #Probably unnecessary, but
  s.suggs[i]    = s.suggs[head]                     #zeroing .sugg zip better.
  s.suggs[head] = CNo(i)

proc allocSz*(s: var Suggestor; head: uint16): Ix {.inline.} =
  result = Ix(s.suggs[head])
  if result == 0:                     #out of space in free list `head`
    s.growArea(head)                  #grow @tail, remap, link in new blocks
    result = Ix(s.suggs[head])
  s.suggs[head] = s.suggs[int(s.suggs[head])]

proc growSuggs*(s: var Suggestor; i: int) =               #"realloc"[i] suggs
  let szIx = te(i).szIx
  let oldI = te(i).sugg
  let curI = s.allocSz(szIx + 1)          #Get slot from next bigger listSz pool
  if int(szIx) > 1:                       #len 0,1 have no lists
    copyMem(addr s.suggs[curI], addr s.suggs[oldI], sizes[szIx] * sizeof(CNo))
  te(i).szIx += 1
  te(i).sugg = curI
  if int(szIx) > 1: s.free(szIx, oldI)    #len 0,1 have no lists

proc addSugg*(s: var Suggestor; d: Word, w: CNo) =        #Add sugg for d
  var i = s.find(d)                                   #Lookup in table, adding
  if i < 0:                                           #key if not present.
    i = addKey(s, i, d, false)
  if int(te(i).len) == 0:                             #len 0,1 are the most..
    te(i).sugg = Ix(w)                                #..common, but these do
    te(i).szIx = 1                                    #..not need any space
  elif int(te(i).len) == 1:                           #..beyond the table since
    let old = CNo(te(i).sugg)                         #..we hijack te(i).sugg
    s.growSuggs(i)                                    #..to point at corpus,
    s.suggs[te(i).sugg] = old                         #..not a list for len 1.
    s.suggs[te(i).sugg+1] = w
  elif int(te(i).len) >= 65534:                       #next inc will wrap
    raise newException(ValueError, "suggestions overflow 16 bits")
  else:
    if int(te(i).len) == sizes[te(i).szIx]:           #Need a longer list
      s.growSuggs(i)                                  #Make room
    s.suggs[te(i).sugg + Ix(te(i).len)] = w           #Add ptr to .meta/.corp
  inc(te(i).len)                                      #register new list length

iterator items*(s: Suggestor, e: TabEnt): CNo =           #Iterate over suggs
  for i in e.sugg ..< e.sugg + e.len:
    yield (if int(e.len) < 2: CNo(e.sugg) else: s.suggs[i])

proc rightSize*(count: Natural): int {.inline.} =
  result = nextPowerOfTwo(count * 8 div 7  +  8)

proc open*(path: string, mode=fmRead, maxDist=3, size=256, refr=""): Suggestor =
  let sz = suggest.rightSize(abs(size))       #Convert size guess to power of 2
  if maxDist < 1: raise newException(ValueError, "Must have maxDist >= 1")
  let tablPath = path & ".tabl"               #Our 5 data files open in about
  let keysPath = path & ".keys"               #10-15 microseconds on Linux.
  let suggPath = path & ".sugg"
  let metaPath = path & ".meta"               #Only need these two for scan.
  let corpPath = path & ".corp"
  result.mode = mode
  let mfop = memfiles.open
  if mode == fmRead:                                #Open existing read only
    if size > 0: result.tabf = mfop(tablPath)
    if size > 0: result.keyf = mfop(keysPath)
    if size > 0: result.sugf = mfop(suggPath)
    result.metf = mfop(metaPath)
    result.corf = mfop(corpPath)
    if refr.len > 0: #HugeTLBfs support uses refr for size; ONLY for fmRead mode
      result.saved = [ result.tabf.size, result.keyf.size, result.sugf.size,
                       result.metf.size, result.corf.size ]
      if size > 0:  result.tabf.size = int(getFileSize(refr & ".tabl"))
      if size > 0:  result.keyf.size = int(getFileSize(refr & ".keys"))
      if size > 0:  result.sugf.size = int(getFileSize(refr & ".sugg"))
      result.metf.size = int(getFileSize(refr & ".meta"))
      result.corf.size = int(getFileSize(refr & ".corp"))
  elif existsFile(tablPath):                        #Open existing Read-Write
    result.tabf = mfop(tablPath, fmReadWrite, allowRemap=true)
    result.keyf = mfop(keysPath, fmReadWrite, allowRemap=true)
    result.sugf = mfop(suggPath, fmReadWrite, allowRemap=true)
    result.metf = mfop(metaPath, fmReadWrite, allowRemap=true)
    result.corf = mfop(corpPath, fmReadWrite, allowRemap=true)
  else:                                             #Create empty files
    result.tabf = mfop(tablPath, fmReadWrite, -1, 0, tBytes(sz), true)
    cast[ptr SuggTab](result.tabf.mem).maxDist = uint8(maxDist)
    result.keyf = mfop(keysPath, fmReadWrite, -1, 0, 1, true)
    result.keyf = mfop(keysPath, fmReadWrite, -1, 0, 1, true)
    result.sugf = mfop(suggPath, fmReadWrite, -1, 0, sizes.len*sizeof(Ix), true)
    result.metf = mfop(metaPath, fmReadWrite, -1, 0, 1, true)
    result.corf = mfop(corpPath, fmReadWrite, -1, 0, 1, true)
  result.table = cast[ptr SuggTab](result.tabf.mem)  #Update Suggestor metadata
  result.tabSz = tSlots(result.tabf.size)            #..inferred from file sizes
  result.suggs = cast[ptr Suggs](result.sugf.mem)
  result.sugSz = Ix(result.sugf.size div sizeof(CNo))
  result.zKeys = if result.keyf.size==1:0 else:result.keyf.size
  result.metas = cast[ptr Metas](result.metf.mem)
  result.nMeta = CNo(if result.metf.size==1: 0
                     else: result.metf.size div sizeof(Meta))
  result.zCorp = if result.corf.size==1:0 else:result.corf.size

proc close*(s: var Suggestor, verb=false, small=false) =
  if verb: echo "size: ", s.tabf.size + int(s.sugSz)*sizeof(CNo) +
                         int(s.zKeys) + int(s.nMeta)*sizeof(Meta) + int(s.zCorp)
  if s.saved[4] != 0:
    s.tabf.size = s.saved[0]; s.keyf.size = s.saved[1]; s.sugf.size = s.saved[2]
    s.metf.size = s.saved[3]; s.corf.size = s.saved[4]
  if not small: memfiles.close(s.tabf)            #This file cannot be trimmed
  if s.mode == fmReadWrite:   #Re-size files so that .size=actually used bytes
    s.keyf.resize(int(s.zKeys))
    s.sugf.resize(int(s.sugSz) * sizeof(CNo))
    s.metf.resize(int(s.nMeta) * sizeof(Meta))
    s.corf.resize(int(s.zCorp))
  if not small: memfiles.close(s.keyf)
  if not small: memfiles.close(s.sugf)
  memfiles.close(s.metf)
  memfiles.close(s.corf)

var diBuf: Word                                   #Global deletes iterator buf
iterator deletes*(w: Word): Word =                #Yield every word w/1 del
  let n = int(w.n)
  if n > 1:                                       #Only an "" delete
    diBuf.n = uint8(n - 1)                        #len is fixed for iteration
    zeroMem(addr diBuf.d, sizeof(diBuf.d))
    copyMem(addr diBuf.d, unsafeAddr w.d[1], n-1) #1st char delete
    yield diBuf
    for i in 1 ..< n - 1:
      copyMem(addr diBuf.d   , unsafeAddr w.d     , i)      #copy from [0..<i]
      copyMem(addr diBuf.d[i], unsafeAddr w.d[i+1], n-1-i)  #copy from [i+1..]
      yield diBuf
    copyMem(addr diBuf.d, unsafeAddr w.d, n - 1)  #Last char delete
    yield diBuf
  elif n == 1:                                    #There is only 1 empty string
    diBuf.n = uint8(0)
    zeroMem(addr diBuf.d, sizeof(diBuf.d))
    yield diBuf

proc hash(w: Word): int {.inline.} = hash(toOpenArray(w.d, 0, int(w.n) - 1))

proc addDeletes*(s: var Suggestor, word: Word, wd: CNo) =
  ## Add strings with up to `maxDist` chars deleted from `word`.
  var adDels   = initSet[Word](1 shl 6)             #Set of all deletions
  var adQueue0 = initSet[Word](4)                   #For depth
  var adQueue1 = adQueue0                           #For depth+1
  var q0 = addr adQueue0
  var q1 = addr adQueue1
  q0[].incl(word)
  for depth in 0 ..< int(s.table[].maxDist) - 1:
    if depth > 0: q1[].clear()                    #Cannot change q0 while itr..
    for w in q0[]:                                #..So, build depth+1 as we go
      for wDel in deletes(w):
        q1[].incl(wDel)                           #Accumulate for next del level
        adDels.incl(wDel)                         #Accumulate for overall dels
    swap q0, q1                                   #Swap ptrs
  for w in q0[]:                                  #Very last level (NOTE ..<)..
    for wDel in deletes(w):                       #..needs no setup for next.
      adDels.incl(wDel)                           #Accumulate for overall dels
  for d in adDels:                                #Add to the delete's sugg seq
    s.addSugg(d, wd)

proc pop[T](s: var HashSet[T]): T {.inline.} = #Rm&Return 1st it-ord elt from s
  for e in s: result = e; break
  s.excl(result)

proc add(s: var Suggestor, word: Word, freq: Count): CNo =
  let zMeta = int(s.nMeta) * sizeof(Meta)
  if zMeta + sizeof(Meta) > s.metf.size:                #Maybe grow .meta file
    s.metf.resize(int(((zMeta + sizeof(Meta)) shl 1)))
    s.metas = cast[ptr Metas](s.metf.mem)
  s.metas[int(s.nMeta)].cnt = freq                      #Init data in .meta
  if int(s.zCorp + Ix(1) + Ix(word.n)) > s.corf.size:   #Maybe grow .corp file
    s.corf.resize(int((s.zCorp + Ix(word.n + 1)) shl 1))
  let eoc = cast[ptr Word](s.corf.mem +% uint(s.zCorp))
  eoc[].n = word.n
  copyMem(addr eoc.d[0], unsafeAddr word.d[0], word.n)  #Init data in .corp
  s.metas[int(s.nMeta)].cix = s.zCorp                   #Update .meta ptr->.corp
  result = s.nMeta                                      #Return corp-word-no
  s.zCorp += Ix(1) + Ix(word.n)                         #Update sizes
  inc(s.nMeta)

proc add*(s: var Suggestor, wrd: string, freq=Count(1)) =
  ## Add a word from some corpus of correct ones with frequencies `freq`.
  let word = toWord(wrd)
  s.table[].total += freq                   #Always update .total
  let i = s.find(word)
  if i >= 0:
    if te(i).corp == 0:                     #Present only as a suggestion
      let cno = s.add(word, freq)
      s.addDeletes(word, cno)               #  Promote to real w/suggs populated
      inc(s.table[].unique)                 #  & register new unique word.
      te(i).corp = 1                        #  with initial frequency `freq`
      te(i).keyI = Ix(cno)
    else:                                   #Already present; just `cnt += freq`
      s.metas[int(te(i).keyI)].cnt += freq
  else:                                     #Absent; common case.
    let j = s.addKey(i, word, true)         #  Add w/no suggs, real word
    let cno = s.add(word, freq)
    s.zKeys -= word.n                       #Trunc .keys
    te(j).corp = 1                          #  with initial frequency `freq`
    te(j).keyI = Ix(cno)
    s.addDeletes(word, cno)                 #  & populate s.tab w/deletes
    inc(s.table[].unique)                   #  & register new unique word.

proc rup(dmaxE: var int; r: Results; d, matches: int) = #Results updater
  var filled = 0                    #Since results are sorted first by distance
  for k in 0..d:                    #and then by other things, we can not add
    filled += r[k].len              #higher distance results and can lower the
    if filled >= matches:           #effective dmax, dmaxE, once we have filled
      dmaxE = d                     #in at least `matches` at a lower [d].
      break                         #Lower dmaxE speeds up future dist calcs.

proc render*(s: Suggestor, rs: Results, matches: int): seq[string] =
  ## render strings sorted by: incr distance,decr frequency,incr discovery order
  var annotated: seq[tuple[d: int, mCnt: Count, cix: Ix]]
  for d in 0 ..< rs.len:            #Collect records for sorting just until we
    for c in rs[d]:                 #..have enough.
      annotated.add( (d, -s.metas[c].cnt, s.metas[c].cix) )
    if annotated.len >= matches:    #NOTE: Could drop both sort() and .meta file
      break                         #..in scan mode if we insist .corp is always
  annotated.sort()                  #..sorted by decreasing frequency.
  result.setLen(min(annotated.len, matches))
  var ss: string
  for i, a in annotated:
    if i == matches: break
    let cw = cast[ucArrCh](s.corf.mem +% a.cix)
    ss.setLen(int(cw[0]))
    copyMem(addr ss[0], addr cw[1], ss.len)
    result[i] = ss

proc suggestions*(s: Suggestor, wrd: string, maxDist: int=3, kind=osa,
                  matches=6): seq[string] =
  ## Return suggested corrections for maybe incorrect word `w`.  `maxDist` may
  ## be usefully less than the `maxDist` used to build Suggestor from a corpus.
  var maxDist = min(maxDist, int(s.table[].maxDist))
  var res: Results; res.setLen(maxDist + 1)
  var added = initSet[CNo](32)
  let p = myersCompile(wrd, 0)                  #For fast distance calc
  let w = toWord(wrd)
  var queue = initSet[Word](32)
  queue.incl(w)                                 #First pop empties this
  while queue.len > 0:
    let q    = queue.pop                        #pop elt q from a set queue
    let entI = s.find(q)
    if entI >= 0:                               #Was present in s.tab
      if te(entI).corp == 1:                    #A word from corpus
        let cno = CNo(te(entI).keyI)
        if cno notin added:
          let cw = cast[ucArrCh](s.corf.mem +% s.metas[cno].cix)
          let d = p.distance(addr cw[1], int(cw[0]), maxDist+1, kind)
          if d <= maxDist:
            res[d].add(cno); added.incl(cno)
            maxDist.rup(res, d, matches)
      for maybe in s.items(te(entI)):           #Check maybe Distance()s
        if maybe in added:                      #Already handled this maybe
          continue
        let cm = cast[ucArrCh](s.corf.mem +% s.metas[maybe].cix)
        let d = p.distance(addr cm[1], int(cm[0]), maxDist+1, kind)
        if d <= maxDist:                        #add close enough maybe's
          res[d].add(maybe); added.incl(maybe)
          maxDist.rup(res, d, matches)
    if int(w.n) - int(q.n) < maxDist:           #Skip if dels too short to pass
      for wDel in deletes(q):                   #incl this q's dels
        queue.incl(wDel)
  result = s.render(res, matches)

proc update*(prefix, input: string; dmax=2, size=32, verbose=false): int =
  ## Build/update input into Suggestor data files stored in `prefix`.* paths.
  let f = memfiles.open(input)
  var s = suggest.open(prefix, fmReadWrite, dmax, size=size)
  var n = 0
  if verbose: echo "Processing input: ", input
  let t0 = epochTime()
  for wordFreq in memSlices(f):
    let cols = ($wordFreq).split()
    if cols[0].len > WHI:                   #Count dropped over-long words
      inc(n)
      continue
    s.add(cols[0], if cols.len > 1: parseFloat(cols[1]) else: 1)
  let dt = epochTime() - t0
  if verbose:
    echo "totalWords: ", s.table[].total, "  unique: ", s.table[].unique
    echo "maxLevD: ", s.table[].maxDist, "  Misspells: ", s.table[].len
    stdout.write "tooLong: ", n, "  dt: ", formatFloat(dt,ffDecimal,3), " sec "
  s.close(verbose)                          #Right-size files when done

proc query*(prefix: string, typos: seq[string], refr="",
            dmax=2, kind=osa, matches=6, verbose=false): int =
  ## Load Suggestor data from `prefix`.* & query suggestions for all `typos`.
  let f00 = tabFinds
  let d00 = totDists
  let t00 = epochTime()
  var s = suggest.open(prefix, refr=refr)   #NOTE: only var so can `.close`
  let dtOp = (epochTime() - t00) * 1e3
  var dtAll = 0.0
  for i in 0 ..< typos.len:
    let f0 = tabFinds
    let d0 = totDists
    let t0 = epochTime()
    let sugg = s.suggestions(typos[i], dmax, kind, matches)
    let dt = (epochTime() - t0) * 1e3
    let df = tabFinds - f0
    let dd = totDists - d0
    dtAll += dt
    stdout.write "  sugg for \"", typos[i], "\""
    if verbose:
      stdout.write " in ", formatFloat(dt, ffDecimal, 4), " ms ",
                   df, " finds ", dd, " dists"
    if sugg.len > 0: stdout.write ":  ", sugg.join(" ")
    echo ""
  if verbose:
    let dp0 = tabFinds - f00
    let dd0 = totDists - d00
    stdout.write formatFloat(dtOp, ffDecimal, 4), " ms to open; ",
                 dp0, " totFind ", dd0, " totDist "
  echo formatFloat(dtAll/float(typos.len), ffDecimal, 4), " ms"
  s.close()

proc suggsScan*(s: Suggestor, typos: seq[string], maxDist: int=3, kind=osa,
                matches=6): seq[seq[string]] =
  let corp = cast[ucArrCh](s.corf.mem)
  var sg = newSeq[Results](typos.len)
  var ps = newSeq[MyersPattern[int]](typos.len)
  var dmaxE = newSeq[int](typos.len)
  for i in 0 ..< typos.len:
    ps[i] = myersCompile(typos[i], 0)
    sg[i].setLen(maxDist + 1)
    dmaxE[i] = maxDist
  var off = 0
  var j = 0
  while off < s.corf.size:
    var n = int(corp[off])                            #First byte is .n
    for i in 0 ..< typos.len:
      let d = ps[i].distance(addr corp[off + 1], n, dmaxE[i] + 1, kind)
      if d < dmaxE[i] + 1:
        sg[i][d].add(CNo(j))
        dmaxE[i].rup(sg[i], d, matches)
    off += n + 1                                      #First byte is .n
    inc(j)
  for sug in sg:
    result.add(s.render(sug, matches))

proc scan*(prefix: string, typos: seq[string], refr="",
           dmax=2, kind=osa, matches=6, verbose=false): int =
  ## Scan prefix.corp (one pass for all typos) to find suggestions.
  let t0 = epochTime()
  var s = suggest.open(prefix, size = -1, refr=refr)  #Only var so can close
  let sg = s.suggsScan(typos, dmax, kind, matches)
  for i, ty in typos:
    stdout.write "  sugg for \"", ty, "\""
    if sg[i].len > 0: stdout.write ":  ", sg[i].join(" ")
    echo ""
  echo formatFloat((epochTime() - t0)*1e3/float(typos.len), ffDecimal, 4), " ms"
  s.close(small=true)

proc makeTypos(path: string, size=6, n=10, deletes=1, outPrefix="typos.") =
  ## Generate ``n` typo files ``outPrefix``* suitable for ``compare`` of size
  ## ``size`` by sampling according to (word,freq) in ``path`` with ``deletes``.
  var words: seq[string]
  var tot = 0.0; var cdf: seq[float]
  for line in system.open(path).lines:
    let cols = line.split
    words.add(cols[0])
    tot += (if cols.len == 2: cols[1].parseFloat else: 1.0); cdf.add(tot)
  for i in 0 ..< n:
    let o = system.open(outPrefix & $i, fmWrite)
    var z = 0
    while z < size:
      var typo = words[cdf.upperBound(rand(cdf[^1]))]
      if typo.len <= deletes: continue    #Need strings > ``deletes`` long
      for d in 0 ..< deletes:
        let k = rand(typo.len - 1)
        typo.delete(k, k)
      o.write typo, "\n"
      inc(z)
    o.close()

proc suggVec*(s: Suggestor, typos: seq[string], dmax=2, kind=osa,
              matches=6): seq[seq[string]] =
  for i in 0 ..< typos.len:
    result.add(s.suggestions(typos[i], dmax, kind, matches))

proc compare*(prefix: string, dir: string, refr="",
              dmax=2, kind=osa, matches=6, verbose=false): int =
  ## A benchmarking call that works with typoGen.py and test-suggest.sh.
  for pkind, path in walkDir(dir):
    var typos: seq[string]
    for line in lines(system.open(path)):
      typos.add(line)
    var t0 = epochTime()
    var s = suggest.open(prefix, size = -1, refr=refr)
    let ss = s.suggsScan(typos, dmax, kind, matches)
    let dtS = (epochTime() - t0) * 1e3
    s.close(small=true)
    t0 = epochTime()
    s = suggest.open(prefix, refr=refr)
    let sq = s.suggVec(typos, dmax, kind, matches)
    let dtQ = (epochTime() - t0) * 1e3
    s.close()
    echo "scan: ", formatFloat(dtS/float(typos.len), ffDecimal, 5),
         " query: ", formatFloat(dtQ/float(typos.len), ffDecimal, 5),
         if ss != sq: " MISMATCH" else: ""
    if verbose:
      for i, lst in ss: echo "  ", typos[i], ": ", lst.join(" ")

proc cpHuge*(paths: seq[string]) =
  ## copy $1 to $2 where $2 is potentially on a Linux hugetlbfs.
  var src = memfiles.open(paths[0])   #FS requires dst.size be a mult of 2M
  let sz2 = (src.size + (1 shl 21) - 1) and not ((1 shl 21) - 1)
  var dst = memfiles.open(paths[1], fmReadWrite, newFileSize=sz2)
  copyMem(dst.mem, src.mem, src.size)
  src.close()
  dst.close()

when defined(Windows):
  proc query2*(): int =
    stderr.write "query2 not implemented for Windows"
    return 1
else:
  import posix    # For TimeVal
  type Rusage* {.importc: "struct rusage", header: "<sys/resource.h>",
                final, pure.} = object
    ru_utime*, ru_stime*: TimeVal
    ru_maxrss*, ru_ixrss*, ru_idrss*, ru_isrss*, ru_minflt*, ru_majflt*, ru_nswap*,
      ru_inblock*, ru_oublock*, ru_msgsnd*, ru_msgrcv*, ru_nsignals*, ru_nvcsw*,
      ru_nivcsw*: clong

  proc getrusage*(who: cint, rusage: ptr Rusage): cint
    {.importc, header: "<sys/resource.h>".}

  proc query2*(prefix: string, typos: seq[string], refr="",
               dmax=2, kind=osa, matches=6): int =
    ## Similar to `query` but per-typo open, close, and measure page faults.
    var s: Suggestor
    var r0, r1: Rusage
    for i in 0 ..< typos.len:
      s = suggest.open(prefix, refr=refr)
      discard getrusage(0, addr r0)
      let sugg = s.suggestions(typos[i], dmax, kind, matches)
      discard getrusage(0, addr r1)
      s.close()
      echo r1.ru_minflt - r0.ru_minflt, " ", typos[i], ": ", sugg.join(" ")

when isMainModule:
  import cligen
  dispatchMulti(
    [ suggest.update, cmdName = "update", help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "input"   : "input corpus file of correct spellings",
      "dmax"    : "max possible distance for queries",
      "size"    : "table size guess",
      "verbose" : "print details about a build/update" } ],
    [ suggest.query, cmdName = "query", help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "dmax"    : "max distance of result from query",
      "kind"    : "kind of distance (osa | lev)",
      "matches" : "max number of matches to print",
      "refr"    : "prefix to refernc files for file sizes",
      "verbose" : "print more details about the query" } ],
    [ suggest.makeTypos, cmdName = "makeTypos" ],
    [ suggest.compare, cmdName = "compare" ],
    [ suggest.cpHuge, cmdName = "cpHuge" ],
    [ suggest.query2, cmdName = "query2" ],
    [ suggest.scan, cmdName = "scan", help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "dmax"    : "max distance of result from query",
      "kind"    : "kind of distance (osa | lev)",
      "matches" : "max number of matches to print",
      "refr"    : "prefix to refernc files for file sizes",
      "verbose" : "print more details about the query" } ])

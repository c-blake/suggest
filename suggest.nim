## A from-scratch re-implementation of Wolf Garbe's Symmetric Delete Algorithm
## for correct spelling suggestions.  It is identical, in essence, to Variation
## 1 of Mor & Fraenkel 1982 "A Hash Code Method for Detecting and Correcting
## Spelling Errors" (DOI 10.1145/358728.358752).
##
## The basic idea is simple: corpus words within edit distance N of a query can
## only be *at most* N units longer or N units shorter and the shorter must be
## derived from deletes of longer strings.  So, build a map from all shortenings
## to all corpus words which generate them.  When queried with a string, we can
## then lookup all possible edits lengthening the query into a corpus word.  We
## also on-the-fly compute all shortenings of a query (& shortenings of
## lengthenings).  From both sets, we filter "maybe within N edits of a corpus
## word" to "actually within N edits".  The filter can be any "distance"
## successfully bounded by N-indels.  This idea is like Monte Carlo numerical
## integration of shapes within easy bounding boxes (but this is deterministic &
## points which pass are reported, not just counted).
##
## While this allows for fast queries, it takes costs much time & space to make
## a lengthenings table.  So, it is a useful strategy if you A) can save the big
## map to disk AND VERY efficiently load it AND/OR B) have MANY queries to
## amortize build costs over.  Rebuilding is lame & "real" DB query latency is
## hostile.  So, this module does an efficient external file format w/5 files to
## mmap&go: A .tabl pointing to (keyIx->varlen.keys,.sugg=varlen[array[CNo]])
## and a .meta(ix,cnt) file pointing to varlen .corp.  varlen[arr[CNo]] is an
## allocation arena with early entries the heads of per-listSz free lists.

import std/[hashes,tables,sets,os,times,memfiles,strutils,algorithm,math,random]
import system/ansi_c
#NOTE: You are NOT intended to understand levenshtein/optimStrAlign here
#      without reading the Hyyro 2003 bit-vector algorithm paper.
type
  MyersPattern*[W] = object       ## Output of myersCompile(var Pattern[W],..)
    m: W                          ## Length of pattern string
    pm: seq[W]                    ## Match vector: pm[ch]&(1<<i)<->A[i]==ch

proc myersCompile*[T,W](pat: openArray[T], dummy: W): MyersPattern[W] =
  ## Compile pattern for either levenshtein or optimStrAlign.
  result.m  = pat.len
  result.pm = newSeq[W](T.high.int)
  for i in 0 ..< pat.len:
    result.pm[pat[i].int] = result.pm[pat[i].int] or (1 shl i)

proc levenshtein*[T,W](p: MyersPattern[W], t: openArray[T], maxDist: W): W =
  ## Myers-Hyyro bit-vec algo (hP..or 1) with early len check against maxDist.
  if p.m == 0: return t.len             #Degenerate case
  let m = p.m
  let n = t.len
  if abs(n - m) >= maxDist:             #Paper tests > & returns +1, but..
    return maxDist                      #..want convention of other procs.
  let mxPBit = 1.W shl (m - 1)          #Bit representing change in last row
  var vP = -1.W                         #Vert change is initially all pos
  var vN = 0.W                          #I.e., no neg.
  var d  = m.W
  for j in 0 ..< n:
    let pmj = p.pm[t[j].int]
    var d0 = (((pmj and vP) + vP) xor vP) or pmj or vN
    var hP = vN or not (d0 or vP)
    let hN = d0 and vP
    if   (hP and mxPBit) != 0: inc d
    elif (hN and mxPBit) != 0: dec d
    hP = (hP shl 1) or 1
    vP = (hN shl 1) or not(d0 or hP)
    vN = d0 and hP
  return min(maxDist, d)

proc optimStrAlign*[T,W](p: MyersPattern[W], t: openArray[T], maxDist: W): W =
  ## Myers-Hyyro bit-vec with early len ck adapted to optimStrAlign. (Hyyro
  ## calls it "Damerau", but it's NOT *unrestricted* transposition Damerau1964
  ## distance.   We call it "osa" to be more clear.)
  if p.m == 0: return t.len             #Degenerate case
  let m = p.m
  let n = t.len
  if abs(n - m) >= maxDist:             #Paper tests > & returns +1, but..
    return maxDist                      #..want convention of other procs.
  let mxPBit = 1.W shl (m - 1)          #Bit representing change in last row
  var vP   = -1.W                       #Vert change is initially all pos
  var vN   = 0.W                        #I.e., no neg.
  var d    = m.W
  var d0   = vN
  var pmj0 = vN
  for j in 0 ..< n:
    let pmj = p.pm[t[j].int]
    d0 = ((((not d0) and pmj) shl 1) and pmj0) or
           (((pmj and vP) + vP) xor vP) or pmj or vN
    pmj0 = pmj
    var hP = vN or not (d0 or vP)
    let hN = d0 and vP
    if   (hP and mxPBit) != 0: inc d
    elif (hN and mxPBit) != 0: dec d
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
  TabEnt* {.packed.} = object     ## Size=13 (should be corp:1 => 12; Nim bug?)
    keyI*: Ix                     ## Possible typo -> suggestion list
    keyN* {.bitsize:  8.}: uint8  ## typo len; (Could be hash w/len in .keys)
    corp* {.bitsize:  1.}: uint8  ## 1-bit flag indicating a real corpus word
    szIx* {.bitsize:  7.}: uint8  ## Index into `sizes` (alloc);
    len*  {.bitsize: 16.}: uint16 ## Length of sugg list (used)
    sugg*: Ix                     ## pointer into file of suggestions lists
  SuggTab* {.packed.} = object    ## 32B header & table
    maxDist*: uint8               ## Within `maxDist` of corpus words
    unique* : CNo                 ## Track unique words add()d
    total*  : float64             ## Track all words add()d
    len*    : Ix                  ## Non-empty entries in table
    tab* : UncheckedArray[TabEnt] ## Hash structured; 12B per entry
  Suggs* = UncheckedArray[CNo]    ## Suggs Arena; 1st 16 are heads of free lists
  Meta* {.packed.} = object
    cix*: Ix                      ## Byte index into .corp for corpus word
    cnt*: Count                   ## Fixed size statistics on the corpus word
  Metas* = UncheckedArray[Meta]   ## Meta data for corpus words
  Suggestor* = object
    cMax*: int                    ## instrumentation: count max probe length
    nFind*: int                   ## instrumentation: count table finds
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
    nMeta*: CNo                   ## USED recs; alloc=metf.size div Meta.sizeof
    corf* : MemFile               ## Backing store for (Count,CorpWord) @CNo's.
    zCorp*: Ix                    ## USED BYTES in corp file; alloc=corf.size
    saved*: array[5, int]         ## HugeTLBfs on Linux rounds st_size.  So, we
                                  ##..get via refr, but munmap wants .saved.
  Results* = seq[seq[CNo]]        ## Per-distance-from-query seq[corpusIx]
  Distance* = enum lev, osa       ## Levenshtein, Optimal String Alignment
  ucArrCh* = ptr UncheckedArray[char] ## Abbrev for raw pointer->char[] cast[]s
let szTabEnt = sizeof(TabEnt)     #Need Nim >=v0.19.9 to have this work

var totDists = 0
proc distance*[W](p: MyersPattern[W], t: ptr char, n: int,
                  maxDist: W, kind = osa): W {.inline.} =
  inc totDists
  if kind == lev: p.levenshtein   toOpenArray(cast[ucArrCh](t), 0, n-1), maxDist
  else          : p.optimStrAlign toOpenArray(cast[ucArrCh](t), 0, n-1), maxDist

proc `+%`[T](p: pointer, i: T): pointer {.inline.} =  #Pointer arithmetic
  cast[pointer](cast[uint](p) + i.uint)

template te[T](i: T): auto = s.table[].tab[i.int]     #TableEnt accessor

#NOTE: .keyI("") != 0 because very first key added is always a real corpus word.
proc empty(e: TabEnt): bool {.inline.} = e.keyI == 0 and e.keyN == 0

proc toWord*(s: string): Word {.inline.} =
  result.n = min(WHI, s.len).uint8
  copyMem result.d[0].addr, s[0].unsafeAddr, result.n

proc `$`*(w: Word): string {.inline.} =
  let n = w.n.int
  result.setLen n
  if n > 0: copyMem result[0].addr, w.d[0].unsafeAddr, n

proc keyData*(s: Suggestor; i: Ix): pointer {.inline.} = #PERSISTENT HASH TABLE
  if te(i).corp == 0:                   #Not corpus word => from .keys
    s.keyf.mem +% te(i).keyI
  else:                                 #corp=>metas->corp; +1 skips .n field
    s.corf.mem +% (s.metas[te(i).keyI.int].cix + 1)

proc toWord*(s: Suggestor; i: Ix): Word {.inline.} =
  result.n = te(i).keyN
  copyMem result.d[0].addr, s.keyData(i), result.n

proc find*(s: var Suggestor, w: Word): int =
  inc s.nFind
  let mask = s.tabSz - 1                #Vanilla linear probe hash search/insert
  var i = w.d.hash and mask             #Initial probe
  let n = w.n.int                       #Length of query word
  var c = 0
  while not te(i).empty:
    inc c
    if i.te.keyN.int == n and
        c_memcmp(s.keyData(i.Ix), w.d[0].unsafeAddr, n.csize_t) == 0:
      return i                          #Len & *only* then bytes equal => Found
    i = (i + 1) and mask                #The linear part of linear probing
  if c > s.cMax:                        #New worst case
    s.cMax = c                          #Print @new max to rate limit messages
    let logM = ln(float(mask + 1))      #See Pittel 1987,"..Probable Largest.."
    let load = s.table[].len.float / float(mask + 1)
    let x = int((logM - 2.5 * ln(logM)) / (load - 1 - ln(load)) + 0.5)
    if c > 2 * x + 14:                  #Lowers chatter for unconcerning cases
      stderr.write "weak hash: ", c, " >> ", x, " loops to find empty (",w,")\n"
  return -i - 1                         #Not Found, return -(insertion point)-1

proc tBytes*(nSlot: int): int {.inline.} = SuggTab.sizeof + nSlot * szTabEnt
proc tSlots*(nByte: int): int {.inline.} = (nByte - SuggTab.sizeof) div szTabEnt

proc growTab*(s: var Suggestor, w: Word): int =
  let oldSz = s.tabSz   #Triple tab[] space, copy old to top, zero new. Might..
  s.tabSz *= 2          #..avoid copy&zero w/more syscalls for marginal speedup.
  s.tabf.resize tBytes(s.tabSz + oldSz)
  s.table = cast[ptr SuggTab](s.tabf.mem)
  copyMem te(s.tabSz).addr, te(0).addr, oldSz*szTabEnt
  zeroMem te(0).addr, s.tabSz * szTabEnt
  for i in s.tabSz ..< s.tabSz + oldSz:             #For each possible OLD slot
    if te(i).empty: continue                        #skip empty
    let J = -s.find(s.toWord i.Ix) - 1              #find spot & copy ent to new
    let j = if J < 0: -J - 1 else: J            #Can falsePos for "" in NEW tab
    copyMem te(j).addr, te(i).addr, szTabEnt    #but cpOver of old non-empty wks
  s.tabf.resize tBytes(s.tabSz)                     #Now shrink file & memory
  return -s.find(w) - 1                   #Return insertion spot in grown table

proc rightSize*(count: Natural): int {.inline.} =
  result = nextPowerOfTwo(count * 8 div 7  +  8)

proc addKey*(s: var Suggestor, mim1: int, w: Word, c: bool): int =
  var i = -mim1 - 1                                 #Add key w/no suggs
  if s.table[].len.int + 1 > 7 * (s.tabSz shr 3):   #> 7/8 => double
    i = s.growTab w
  s.table[].len += 1
  if int(s.zKeys + w.n.Ix) > s.keyf.size:           #Maybe grow .keys file
    s.keyf.resize int(s.zKeys + w.n.Ix) shl 1
  copyMem s.keyf.mem +% s.zKeys, w.d[0].unsafeAddr, w.n
  te(i).keyI = s.zKeys                              #Point table at new key
  te(i).keyN = w.n
  s.zKeys += w.n.Ix
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
  if s.sugSz.int * CNo.sizeof > s.sugf.size:    #Need more file space:
    s.sugf.resize (s.sugSz.int * CNo.sizeof) shl 1  #Double, remap and fix..
    s.suggs = cast[ptr Suggs](s.sugf.mem)           #..convenience field.
  s.suggs[head] = oldSz.CNo                     #head = oldSz & thread rest
  if grows[head] >= 2:                          #New space from OS is all 0
    var ix = oldSz                              #NOTE: This free-list threading
    while ix < oldSz + Ix(grows[head] - 2):     #..is inactive with grows[]=1.
      s.suggs[ix] = CNo(ix + sizes[head].Ix)
      ix += sizes[head].Ix

proc free*(s: var Suggestor; head: uint16, i: Ix) {.inline.} =
  zeroMem s.suggs[i].addr, Ix.sizeof*sizes[head]     #Probably unnecessary, but
  s.suggs[i]    = s.suggs[head]                      #zeroing .sugg zips better.
  s.suggs[head] = i.CNo

proc allocSz*(s: var Suggestor; head: uint16): Ix {.inline.} =
  result = s.suggs[head].Ix
  if result == 0:                     #out of space in free list `head`
    s.growArea head                   #grow @tail, remap, link in new blocks
    result = s.suggs[head].Ix
  s.suggs[head] = s.suggs[s.suggs[head].int]

proc growSuggs*(s: var Suggestor; i: int) =               #"realloc"[i] suggs
  let szIx = te(i).szIx
  let oldI = te(i).sugg
  let curI = s.allocSz(szIx + 1)          #Get slot from next bigger listSz pool
  if szIx.int > 1:                        #len 0,1 have no lists
    copyMem s.suggs[curI].addr, s.suggs[oldI].addr, sizes[szIx]*CNo.sizeof
  te(i).szIx += 1
  te(i).sugg = curI
  if szIx.int > 1: s.free szIx, oldI      #len 0,1 have no lists

proc addSugg*(s: var Suggestor; d: Word, w: CNo) =        #Add sugg for d
  var i = s.find d                                    #Lookup in table, adding
  if i < 0:                                           #key if not present.
    i = s.addKey(i, d, false)
  if te(i).len.int == 0:                              #len 0,1 are the most..
    te(i).sugg = w.Ix                                 #..common, but these do
    te(i).szIx = 1                                    #..not need any space
  elif te(i).len.int == 1:                            #..beyond the table since
    let old = te(i).sugg.CNo                          #..we hijack te(i).sugg
    s.growSuggs i                                     #..to point at corpus,
    s.suggs[te(i).sugg] = old                         #..not a list for len 1.
    s.suggs[te(i).sugg + 1] = w
  elif te(i).len.int >= 65534:                        #next inc will wrap
    raise newException(ValueError, "suggestions overflow 16 bits")
  else:
    if te(i).len.int == sizes[te(i).szIx]:            #Need a longer list
      s.growSuggs i                                   #Make room
    s.suggs[te(i).sugg + te(i).len.Ix] = w            #Add ptr to .meta/.corp
  inc te(i).len                                       #register new list length

iterator items*(s: Suggestor, e: TabEnt): CNo =           #Iterate over suggs
  for i in e.sugg ..< e.sugg + e.len:
    yield (if e.len.int < 2: e.sugg.CNo else: s.suggs[i])

proc open*(path: string, mode=fmRead, maxDist=3, size=256, refr=""): Suggestor =
  let sz = suggest.rightSize(abs size)        #Convert size guess to power of 2
  if maxDist < 1: raise newException(ValueError, "Must have maxDist >= 1")
  let tablPath = path & ".tabl"               #Our 5 data files open in about
  let keysPath = path & ".keys"               #10-15 microseconds on Linux.
  let suggPath = path & ".sugg"
  let metaPath = path & ".meta"               #Only need these two for scan.
  let corpPath = path & ".corp"
  result.mode = mode
  let mfop = memfiles.open
  if mode == fmRead:                                #Open existing read only
    if size > 0: result.tabf = mfop tablPath
    if size > 0: result.keyf = mfop keysPath
    if size > 0: result.sugf = mfop suggPath
    result.metf = mfop metaPath
    result.corf = mfop corpPath
    if refr.len > 0: #HugeTLBfs support uses refr for size; ONLY for fmRead mode
      result.saved = [ result.tabf.size, result.keyf.size, result.sugf.size,
                       result.metf.size, result.corf.size ]
      if size > 0: result.tabf.size = getFileSize(refr & ".tabl").int
      if size > 0: result.keyf.size = getFileSize(refr & ".keys").int
      if size > 0: result.sugf.size = getFileSize(refr & ".sugg").int
      result.metf.size = getFileSize(refr & ".meta").int
      result.corf.size = getFileSize(refr & ".corp").int
  elif fileExists tablPath:                         #Open existing Read-Write
    result.tabf = mfop(tablPath, fmReadWrite, allowRemap=true)
    result.keyf = mfop(keysPath, fmReadWrite, allowRemap=true)
    result.sugf = mfop(suggPath, fmReadWrite, allowRemap=true)
    result.metf = mfop(metaPath, fmReadWrite, allowRemap=true)
    result.corf = mfop(corpPath, fmReadWrite, allowRemap=true)
  else:                                             #Create empty files
    result.tabf = mfop(tablPath, fmReadWrite, -1, 0, tBytes sz, true)
    cast[ptr SuggTab](result.tabf.mem).maxDist = maxDist.uint8
    result.keyf = mfop(keysPath, fmReadWrite, -1, 0, 1, true)
    result.keyf = mfop(keysPath, fmReadWrite, -1, 0, 1, true)
    result.sugf = mfop(suggPath, fmReadWrite, -1, 0, sizes.len*Ix.sizeof, true)
    result.metf = mfop(metaPath, fmReadWrite, -1, 0, 1, true)
    result.corf = mfop(corpPath, fmReadWrite, -1, 0, 1, true)
  result.table = cast[ptr SuggTab](result.tabf.mem)  #Update Suggestor metadata
  result.tabSz = tSlots result.tabf.size             #..inferred from file sizes
  result.suggs = cast[ptr Suggs](result.sugf.mem)
  result.sugSz = Ix(result.sugf.size div CNo.sizeof)
  result.zKeys = if result.keyf.size==1: 0 else: result.keyf.size
  result.metas = cast[ptr Metas](result.metf.mem)
  result.nMeta = CNo(if result.metf.size==1: 0
                     else: result.metf.size div Meta.sizeof)
  result.zCorp = if result.corf.size==1: 0 else: result.corf.size

proc close*(s: var Suggestor, verb=false, small=false) =
  if verb: echo "size: ", s.tabf.size + s.sugSz.int * CNo.sizeof +
                          s.zKeys.int + s.nMeta.int * Meta.sizeof + s.zCorp.int
  if s.saved[4] != 0:
    s.tabf.size = s.saved[0]; s.keyf.size = s.saved[1]; s.sugf.size = s.saved[2]
    s.metf.size = s.saved[3]; s.corf.size = s.saved[4]
  if not small: memfiles.close s.tabf             #This file cannot be trimmed
  if s.mode == fmReadWrite:   #Re-size files so that .size=actually used bytes
    s.keyf.resize s.zKeys.int
    s.sugf.resize s.sugSz.int * CNo.sizeof
    s.metf.resize s.nMeta.int * Meta.sizeof
    s.corf.resize s.zCorp.int
  if not small: memfiles.close s.keyf
  if not small: memfiles.close s.sugf
  memfiles.close s.metf
  memfiles.close s.corf

var diBuf: Word                                   #Global deletes iterator buf
iterator deletes*(w: Word): Word =                #Yield every word w/1 del
  let n = w.n.int
  if n > 1:                                       #Only an "" delete
    diBuf.n = uint8(n - 1)                        #len is fixed for iteration
    zeroMem diBuf.d.addr, diBuf.d.sizeof
    copyMem diBuf.d.addr, w.d[1].unsafeAddr, n-1  #1st char delete
    yield diBuf
    for i in 1 ..< n - 1:
      copyMem diBuf.d   .addr, w.d.unsafeAddr     , i       #copy from [0..<i]
      copyMem diBuf.d[i].addr, w.d[i+1].unsafeAddr, n-1-i   #copy from [i+1..]
      yield diBuf
    copyMem diBuf.d.addr, w.d.unsafeAddr, n-1     #Last char delete
    yield diBuf
  elif n == 1:                                    #There is only 1 empty string
    diBuf.n = 0.uint8
    zeroMem diBuf.d.addr, diBuf.d.sizeof
    yield diBuf

proc hash(w: Word): int {.inline.} = toOpenArray(w.d, 0, w.n.int - 1).hash

proc addDeletes*(s: var Suggestor, word: Word, wd: CNo) =
  ## Add strings with up to `maxDist` chars deleted from `word`.
  var adDels   = initHashSet[Word](1 shl 6)         #Set of all deletions
  var adQueue0 = initHashSet[Word](4)               #For depth
  var adQueue1 = adQueue0                           #For depth+1
  var q0 = adQueue0.addr
  var q1 = adQueue1.addr
  q0[].incl word
  for depth in 0 ..< s.table[].maxDist.int - 1:
    if depth > 0: q1[].clear                      #Cannot change q0 while itr..
    for w in q0[]:                                #..So, build depth+1 as we go
      for wDel in w.deletes:
        q1[].incl wDel                            #Accumulate for next del level
        adDels.incl wDel                          #Accumulate for overall dels
    swap q0, q1                                   #Swap ptrs
  for w in q0[]:                                  #Very last level (NOTE ..<)..
    for wDel in w.deletes:                        #..needs no setup for next.
      adDels.incl wDel                            #Accumulate for overall dels
  for d in adDels:                                #Add to the delete's sugg seq
    s.addSugg d, wd

proc add(s: var Suggestor, word: Word, freq: Count): CNo =
  let zMeta = s.nMeta.int * Meta.sizeof
  if zMeta + Meta.sizeof > s.metf.size:                 #Maybe grow .meta file
    s.metf.resize (zMeta + Meta.sizeof) shl 1
    s.metas = cast[ptr Metas](s.metf.mem)
  s.metas[s.nMeta.int].cnt = freq                       #Init data in .meta
  let zCorpNew = int(s.zCorp + Ix(word.n + 1))
  if zCorpNew > s.corf.size:                            #Maybe grow .corp file
    s.corf.resize zCorpNew shl 1
  let eoc = cast[ptr Word](s.corf.mem +% s.zCorp)       #End Of Corp as a Word
  eoc[].n = word.n
  copyMem eoc.d[0].addr, word.d[0].unsafeAddr, word.n   #Init data in .corp
  s.metas[s.nMeta.int].cix = s.zCorp                    #Update .meta ptr->.corp
  result = s.nMeta                                      #Return corp-word-no
  s.zCorp = zCorpNew.Ix                                 #Update sizes
  inc s.nMeta

proc add*(s: var Suggestor, wrd: string, freq=Count(1)) =
  ## Add a word from some corpus of correct ones with frequencies `freq`.
  let word = wrd.toWord
  s.table[].total += freq                   #Always update .total
  let i = s.find word
  if i >= 0:
    if te(i).corp == 0:                     #Present only as a suggestion
      let cno = s.add(word, freq)
      s.addDeletes(word, cno)               #  Promote to real w/suggs populated
      inc s.table[].unique                  #  & register new unique word.
      te(i).corp = 1                        #  with initial frequency `freq`
      te(i).keyI = cno.Ix
    else:                                   #Already present; just `cnt += freq`
      s.metas[te(i).keyI.int].cnt += freq
  else:                                     #Absent; common case.
    let j = s.addKey(i, word, true)         #  Add w/no suggs, real word
    let cno = s.add(word, freq)
    s.zKeys -= word.n                       #Trunc .keys
    te(j).corp = 1                          #  with initial frequency `freq`
    te(j).keyI = cno.Ix
    s.addDeletes word, cno                  #  & populate s.tab w/deletes
    inc s.table[].unique                    #  & register new unique word.

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
  annotated.sort                    #..sorted by decreasing frequency.
  result.setLen min(annotated.len, matches)
  var ss: string
  for i, a in annotated:
    if i == matches: break
    let cw = cast[ucArrCh](s.corf.mem +% a.cix)
    ss.setLen cw[0].int
    copyMem ss[0].addr, cw[1].addr, ss.len
    result[i] = ss

proc suggestions*(s: var Suggestor, wrd: string, maxDist: int=3, kind=osa,
                  matches=6): seq[string] =
  ## Return suggested corrections for maybe incorrect word `w`.  `maxDist` may
  ## be usefully less than the `maxDist` used to build Suggestor from a corpus.
  var maxDist = min(maxDist, s.table[].maxDist.int)
  var res: Results; res.setLen(maxDist + 1)
  var added = initHashSet[CNo](32)
  let p = myersCompile(wrd, 0)                  #For fast distance calc
  let w = wrd.toWord
  var queue = initHashSet[Word](32)
  queue.incl w                                  #First pop empties this
  while queue.len > 0:
    let q    = queue.pop                        #pop elt q from a set queue
    let entI = s.find q
    if entI >= 0:                               #Was present in s.tab
      if te(entI).corp == 1:                    #A word from corpus
        let cno = te(entI).keyI.CNo
        if cno notin added:
          let cw = cast[ucArrCh](s.corf.mem +% s.metas[cno].cix)
          let d = p.distance(cw[1].addr, cw[0].int, maxDist+1, kind)
          if d <= maxDist:
            res[d].add cno; added.incl cno
            maxDist.rup res, d, matches
      for maybe in s.items te(entI):            #Check maybe Distance()s
        if maybe in added:                      #Already handled this maybe
          continue
        let cm = cast[ucArrCh](s.corf.mem +% s.metas[maybe].cix)
        let d = p.distance(cm[1].addr, cm[0].int, maxDist+1, kind)
        if d <= maxDist:                        #add close enough maybe's
          res[d].add maybe; added.incl maybe
          maxDist.rup res, d, matches
    if w.n.int - q.n.int < maxDist:             #Skip if dels too short to pass
      for wDel in q.deletes:                    #incl this q's dels
        queue.incl wDel
  result = s.render(res, matches)

proc update*(prefix, input: string; dmax=2, size=32, dlm=' ',verbose=false):int=
  ## Build/update input into Suggestor data files in `prefix`.\* paths.
  let f = memfiles.open input
  var s = suggest.open(prefix, fmReadWrite, dmax, size=size)
  var n = 0
  if verbose: echo "Processing input: ", input
  let t0 = epochTime()
  for wordFreq in f.memSlices:
    let cols = split($wordFreq, dlm)
    if cols[0].len > WHI:                   #Count dropped over-long words
      inc n
      continue
    s.add(cols[0], if cols.len > 1: parseFloat(cols[1]) else: 1)
  let dt = epochTime() - t0
  if verbose:
    echo "totalWords: ", s.table[].total, "  unique: ", s.table[].unique
    echo "maxLevD: ", s.table[].maxDist, "  Misspells: ", s.table[].len
    stdout.write "tooLong: ", n, "  dt: ", formatFloat(dt,ffDecimal,3), " sec "
  s.close verbose                           #Right-size files when done

proc query*(prefix: string, typos: seq[string], refr="",
            dmax=2, kind=osa, matches=6, verbose=false, n=1): int =
  ## Load Suggestor data from `prefix`.\* & query suggestions for `typos`
  var dp0 = 0
  var dd0 = 0
  let t00 = epochTime()
  var s = suggest.open(prefix, refr=refr)   #NOTE: only var so can `.close`
  let dtOp = (epochTime() - t00) * 1e3
  var dtAll = 0.0
  for i in 0 ..< typos.len:
    let f0 = s.nFind
    let d0 = totDists
    let t0 = epochTime()
    for reps in 1..n-1: discard s.suggestions(typos[i], dmax, kind, matches)
    let sugg = s.suggestions(typos[i], dmax, kind, matches)
    let dt = (epochTime() - t0) * 1e3
    let df = s.nFind - f0 ; dp0 += df
    let dd = totDists - d0; dd0 += dd
    dtAll += dt
    stdout.write "  sugg for \"", typos[i], "\""
    if verbose:
      stdout.write " in ", formatFloat(dt, ffDecimal, 4), " ms ",
                   df, " s.nFind ", dd, " dists"
    if sugg.len > 0: stdout.write ":  ", sugg.join(" ")
    echo ""
  if verbose:
    stdout.write formatFloat(dtOp, ffDecimal, 4), " ms to open; ",
                 dp0, " totFind ", dd0, " totDist "
  echo formatFloat(dtAll/typos.len.float, ffDecimal, 4), " ms"
  s.close

proc suggsScan*(s: Suggestor, typos: seq[string], maxDist: int=3, kind=osa,
                matches=6): seq[seq[string]] =
  let corp = cast[ucArrCh](s.corf.mem)
  var sg = newSeq[Results](typos.len)
  var ps = newSeq[MyersPattern[int]](typos.len)
  var dmaxE = newSeq[int](typos.len)
  for i in 0 ..< typos.len:
    ps[i] = myersCompile(typos[i], 0)
    sg[i].setLen maxDist + 1
    dmaxE[i] = maxDist
  var off = 0
  var j = 0
  while off < s.corf.size:
    var n = corp[off].int                             #First byte is .n
    for i in 0 ..< typos.len:
      let d = ps[i].distance(corp[off + 1].addr, n, dmaxE[i] + 1, kind)
      if d < dmaxE[i] + 1:
        sg[i][d].add j.CNo
        dmaxE[i].rup sg[i], d, matches
    off += n + 1                                      #First byte is .n
    inc j
  for sug in sg:
    result.add s.render(sug, matches)

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
  echo formatFloat((epochTime() - t0)*1e3/typos.len.float, ffDecimal, 4), " ms"
  s.close(small=true)

proc makeTypos(path: string, size=6, n=10, deletes=1, outPrefix="typos.") =
  ## Make ``n`` typo files w/``size`` entries in ``outPrefix``* as input for
  ## ``compare``. Sample via w,frq in ``path`` &do ``deletes`` dels.
  var words: seq[string]
  var tot = 0.0; var cdf: seq[float]
  for line in system.open(path).lines:
    let cols = line.split
    words.add cols[0]
    tot += (if cols.len == 2: cols[1].parseFloat else: 1.0)
    cdf.add tot
  for i in 0 ..< n:
    let o = system.open(outPrefix & $i, fmWrite)
    var z = 0
    while z < size:
      var typo = words.sample cdf
      if typo.len <= deletes: continue    #Need strings > ``deletes`` long
      for d in 0 ..< deletes:
        let k = rand(typo.len - 1)
        typo.delete(k, k)
      o.write typo, "\n"
      inc z
    o.close

proc suggVec*(s: var Suggestor, typos: seq[string], dmax=2, kind=osa,
              matches=6): seq[seq[string]] =
  for i in 0 ..< typos.len:
    result.add s.suggestions(typos[i], dmax, kind, matches)

proc compare*(prefix: string, dir: string, refr="",
              dmax=2, kind=osa, matches=6, verbose=false): int =
  ## A benchmarking call that works with makeTypos and test-suggest.sh.
  for pkind, path in dir.walkDir:
    var typos: seq[string]
    for line in system.open(path).lines:
      typos.add line
    var t0 = epochTime()
    var s = suggest.open(prefix, size = -1, refr=refr)
    let ss = s.suggsScan(typos, dmax, kind, matches)
    let dtS = (epochTime() - t0) * 1e3
    s.close(small=true)
    t0 = epochTime()
    s = suggest.open(prefix, refr=refr)
    let sq = s.suggVec(typos, dmax, kind, matches)
    let dtQ = (epochTime() - t0) * 1e3
    s.close
    echo "scan: ", formatFloat(dtS/typos.len.float, ffDecimal, 5),
         " query: ", formatFloat(dtQ/typos.len.float, ffDecimal, 5),
         if ss != sq: " MISMATCH" else: ""
    if verbose:
      for i, lst in ss: echo "  ", typos[i], ": ", lst.join(" ")

proc cpHuge*(paths: seq[string]) =
  ## cpy paths[0] to paths[1] when paths[1] may be on a Linux hugetlbfs.
  var src = memfiles.open paths[0]
  let rnd = (1 shl 21) - 1        #FS requires dst.size be a mult of 2M
  var dst = memfiles.open(paths[1], fmReadWrite,
                          newFileSize = (src.size + rnd) and not rnd)
  copyMem dst.mem, src.mem, src.size
  src.close
  dst.close

when defined Windows:
  proc iquery*(): int =
    stderr.write "iquery not implemented for Windows"
    return 1  #PR welcome someone wants to impl/test GetProcessMemoryInfo
else:
  import posix
  proc iquery*(prefix: string, typos: seq[string], refr="",
               dmax=2, kind=osa, matches=6): int =
    ## Like `query` but per-typo open, close, and page fault instrumented.
    var s: Suggestor
    var r0, r1: Rusage
    for i in 0 ..< typos.len:
      s = suggest.open(prefix, refr=refr)
      getrusage(RUSAGE_SELF, r0.addr)
      let sugg = s.suggestions(typos[i], dmax, kind, matches)
      getrusage(RUSAGE_SELF, r1.addr)
      s.close
      echo r1.ru_minflt - r0.ru_minflt, " ", typos[i], ": ", sugg.join(" ")

when isMainModule:
  import cligen
  include cligen/mergeCfgEnv      #Needs very recent cligen
  dispatchMulti([ "multi", doc = "" ],
    [ suggest.update, help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "input"   : "input corpus file of correct spellings",
      "dmax"    : "max possible distance for queries",
      "size"    : "table size guess",
      "verbose" : "print details about a build/update" } ],
    [ suggest.query, help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "dmax"    : "max distance of result from query",
      "kind"    : "kind of distance (osa | lev)",
      "matches" : "max number of matches to print",
      "refr"    : "prefix to refernc files for file sizes",
      "verbose" : "print more details about the query" } ],
    [ suggest.scan, help = {
      "prefix"  : "path prefix for .keys, etc. data files",
      "dmax"    : "max distance of result from query",
      "kind"    : "kind of distance (osa | lev)",
      "matches" : "max number of matches to print",
      "refr"    : "prefix to refernc files for file sizes",
      "verbose" : "print more details about the query" } ],
    [ suggest.compare ], [ makeTypos ], [ iquery ], # A few for benchmarking
    [ cpHuge, usage = "$command $args\n${doc}" ])

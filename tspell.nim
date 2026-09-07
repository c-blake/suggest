## Standalone Nim for Peter Schneider-Kamp, "Approximate Dictionary Searching at
## a Scale using Ternary Search Trees and Implicit Levenshtein Automata", ICSOFT
## 2022, pp.657-662.  Fixes 1 NASTY bug, 2 obvious bugs & 2 major omissions from
## the paper.  This version also works off a [0]=nil file-friendly Node pool.
import std/[tables, sugar, algorithm]
type Node* = object     ## 20B TernST Node: no GC hdr, no malloc/node.
  c*,d,e,f: char ## Byte,Pads; 29bit link ptrs could do 16B;Slower&smaller alph
  v*: float32   ## Here used for freq; Called 'v' since can be any user value!=0
  l*, m*, r*: uint32 ## 3 links in ternary search tree
type Nodes* = seq[Node] ## Pool ~10% faster; Main cost=bad specul/branch.

proc add*(ns: var Nodes; n:uint32; s:cstring; len,o: int; v:float32): uint32 =
  if n == 0:  # [n].c == => `<`, `>` surely fail => Skip & inline rest here
    let n = ns.len.uint32; ns.add Node(c: s[o])     # How we alloc a Node
    if o+1 < len: ns[n].m = ns.add(0, s, len, o+1, v); result = n
    else        : ns[n].v = v                        ; result = n
  elif s[o] < ns[n].c: ns[n].l = ns.add(ns[n].l, s, len, o  , v); result = n
  elif s[o] > ns[n].c: ns[n].r = ns.add(ns[n].r, s, len, o  , v); result = n
  elif o+1 < len     : ns[n].m = ns.add(ns[n].m, s, len, o+1, v); result = n
  else: ns[n].v = v; result = n       # Paper bug in Above line; =~ n.r.add

type Res* = Table[string, (int, float32)] # Results; TODO? Small `matches`=>seq

proc add(r: var Res, e: string, t: int, v: float32, dMx: int) =
  let d = dMx - t                               # Paper omit: Dist unreported
  try            : r[e][0] = r[e][0].min d      # Paper omit: Many dupes without
  except KeyError: r[e]    = (d, v)             #..some lookup filter like this.

proc near*(r: var Res, ns: Nodes, n: uint32, s: string, t: int; v=0f32, o=0,
           hasV=false, w=false, d=false, e: var string, dMx=t) =
  ## Algo3: GET(n,s,t,v,w,d,e). w=(w)ildcard flag. d=edit explore already (d)one
  ## for this part of path. Adapted to use Node pool, allocate only once/hit by
  ## updating `o` (offset) into `s` (not changing it) & by `var e` updates.
  if o == s.len and not w and hasV: r.add e, t, v, dMx
  if n != 0 and (o < s.len or w):
    let b = ns[n]       # Current Node (B)ody/o(B)j
    if w or s[o] < b.c: r.near ns,b.l, s, t, 0, o, false, w, true, e, dMx
    if w or s[o] > b.c: r.near ns,b.r, s, t, 0, o, false, w, true, e, dMx
    if w or s[o] == b.c:
      e.add b.c
      r.near ns,b.m, s, t, ns[n].v, o+int(not w), b.v != 0, false, false, e, dMx
      e.setLen e.len - 1
  if not d and t >= 1 and not w: # May edit; Paper bug `¬c` fixed to `¬d`
    if n != 0   : r.near ns,n, s, t-1, 0, o  , false, true, false, e, dMx # Ins
    if o < s.len: r.near ns,n, s, t-1, 0, o+1, hasV, false, false, e, dMx # Del
    if n != 0 and o < s.len:                # ^^^^- Always false=NASTY Paper bug
                  r.near ns,n, s, t-1, 0, o+1, false, true, false, e, dMx # Sub

proc ord*(r:Res):auto = result=collect(for w,v in r: (v[0],-v[1],w));result.sort
  ## Sort first by increasing dist THEN decreasing weight THEN increasing alpha.

when isMainModule:
  import std/[strutils, times, syncio], cligen, cligen/[mslice, osUt]
  proc memchr(s:cstring, c:char, n:int): pointer {.importc, header:"string.h".}
  proc tspell(typos: seq[string], freqs: string, z=4, dMx=2, matches=5, verb=2)=
    ## `suggest`-like spell-check. `freqs` format: Word<SingleSpace>IntCount\\n.
    proc `$`(r: Res): string =
      var ws: seq[string]
      for (_,_,w) in r.ord: (if ws.len == matches: break else: ws.add w)
      ws.join " "
    var ns = newSeqOfCap[Node](z)     # ix 0 reserved=nil sentinel; Real[] >=1
    var t = 0u32                      # root ix (0 = empty tree)
    let t0 = epochTime() # Building is 20ms affair; Could save via mmap-alloc
    for (cs, n) in freqs.getDelims:   # Simple input format w/exactly 1-space
      let p = memchr(cs, ' ', n); if p.isNil: quit "Non 2-col fmt `freqs`", 1
      let m = p -! cs
      t = ns.add(t, cs, m,0, MSlice(mem: p+!1, len: n-m-1).parseFloat.float32)
    let t1 = epochTime()
    var e = newString(64); var r: Res # Re-use hot memory from one typo to next
    for typo in typos:
      r.clear # To allow mem re-use, `r` & `e` reset a part of call convention.
      for d in 1..dMx:
        e.setLen 0; r.near(ns,t, typo, d, e=e)
        if matches > 0 and r.len >= matches: break
      if verb > 0: echo typo," (",r.len,"): ",(if verb > 1: $r else: "")
    let t2 = epochTime(); template ff3(v): untyped = formatFloat(v,ffDecimal,3)
    if verb==0: stderr.write "nodes: ",ns.len-1," bytes: ",ns.len*Node.sizeof,
                      " build/file: ",ff3(1e3*(t1 - t0)                )," ms",
                      " query/typo: ",ff3(1e3*(t2 - t1)/typos.len.float)," ms\n"
  include cligen/mergeCfgEnv
  dispatch tspell, help={ "typos": "list of words for which to gen suggestions",
    "freqs": "path to integer-weighted dictionary",
    "z": "pre-siZe this many nodes", "dMx": "max allowed edit distance",
    "matches": "max suggestion count",
    "verb": "0:Just timing; 1:Just counts; 2:Full report"}

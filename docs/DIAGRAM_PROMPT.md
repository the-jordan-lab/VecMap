# Diagram Generation Prompt

For creating Figure 2: NumPy vectorization vs naive Python loops

## Prompt for diagram tools (BioRender, draw.io, Excalidraw, etc.)

**Title:** "Exact-match alignment in NumPy vs. naïve Python loops"

**Scene layout:** Two horizontal panels sitting side-by-side.

### Left panel – "VecMap NumPy vectorization"

* Show a read sequence (rectangle "READ") and *N* candidate reference windows (stacked rectangles "REF slice 1…N") entering a *single* NumPy broadcast/comparison block.
* Inside the block, depict a 2-D Boolean matrix being generated (mismatch matrix).
* Immediately downstream, show a `sum(axis=1)` reduction producing a 1-D vector of mismatch counts, followed by an `argmin()` selector pointing to the best hit.
* Annotate: "all candidates scored in **one** array operation; work done in C under the hood".

### Right panel – "Naïve Python loop"

* Show the same read rectangle feeding into a **for-loop** symbol iterating over "candidate i".
* For each iteration, a small mismatch-counting loop icon (nested loop) compares bases one-by-one.
* Place a red clock icon beside the outer loop to stress cumulative interpreter overhead.

### Bridging note

* Draw a dotted line between the two panels labelled "Same algorithmic logic – different implementation strategy."
* Add a discreet footnote at the bottom: "Python-level loops introduce interpreter overhead; vectorization moves the heavy lifting to highly-optimized C code – experts will recognise why this matters."

**Colour palette:** neutral greys for boxes, soft blue highlight for NumPy array, light orange for Python loops (avoid bright, distracting colours).

**Output size:** 1600 × 800 px SVG preferred (any vector format acceptable). 
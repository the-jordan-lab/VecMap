# VecMap bioRxiv Submission Package

This directory contains the LaTeX source files for the VecMap manuscript prepared for bioRxiv submission.

## Files Included

- `manuscript_vecmap_biorxiv.tex` - Main LaTeX manuscript file
- `references.bib` - Bibliography file with all citations
- `docs/figures/` - Directory containing all figure files:
  - `figure1_performance_comparison.pdf` - Performance benchmarks
  - `figure2_vectorization_impact.pdf` - Vectorization speedup analysis
  - `figure3_crispr_detail.pdf` - CRISPR guide detection performance
  - Additional supplementary figures

## Compilation Instructions

1. **Using pdflatex + bibtex (traditional method):**
   ```bash
   pdflatex manuscript_vecmap_biorxiv
   bibtex manuscript_vecmap_biorxiv
   pdflatex manuscript_vecmap_biorxiv
   pdflatex manuscript_vecmap_biorxiv
   ```

2. **Using latexmk (recommended):**
   ```bash
   latexmk -pdf manuscript_vecmap_biorxiv.tex
   ```

3. **Using Overleaf:**
   - Upload all files maintaining the directory structure
   - Compile with pdfLaTeX

## bioRxiv Submission Checklist

- [x] Line numbers enabled (using `\linenumbers`)
- [x] Double spacing (1.5x via `\renewcommand{\baselinestretch}{1.5}`)
- [x] Figures included as PDF files
- [x] Complete author information and affiliations
- [x] Abstract under 350 words
- [x] Keywords provided
- [x] Data availability statement included
- [x] Acknowledgments section
- [x] References in standard format

## Key Highlights

- **Performance**: 42,027 reads/second in pure Python
- **CRISPR Applications**: >1M reads/second for guide detection
- **Memory Efficient**: Only 22MB for transcriptome alignment
- **Educational Value**: Simple 88-line implementation

## Abstract Word Count

Current abstract: ~200 words (well under bioRxiv's 350-word limit)

## Corresponding Author

James M. Jordan  
Department of Biological Science  
Florida State University  
Email: jmjordan@fsu.edu, jmjordan@bio.fsu.edu, or jim@jordanlab.org

## Repository

Code available at: https://github.com/the-jordan-lab/VecMap

## License

The manuscript is available under CC-BY 4.0 license for bioRxiv. 
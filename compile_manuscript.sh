#!/bin/bash

# Compile VecMap manuscript for bioRxiv

echo "Compiling VecMap manuscript..."

# Check if pdflatex is available
if ! command -v pdflatex &> /dev/null; then
    echo "Error: pdflatex not found. Please install a LaTeX distribution."
    exit 1
fi

# Compile the document
pdflatex -interaction=nonstopmode manuscript_vecmap_biorxiv.tex
if [ $? -eq 0 ]; then
    bibtex manuscript_vecmap_biorxiv
    pdflatex -interaction=nonstopmode manuscript_vecmap_biorxiv.tex
    pdflatex -interaction=nonstopmode manuscript_vecmap_biorxiv.tex
    echo "✅ Compilation successful! Output: manuscript_vecmap_biorxiv.pdf"
else
    echo "❌ Compilation failed. Check the log file for errors."
    exit 1
fi

# Clean up auxiliary files (optional)
# rm -f *.aux *.bbl *.blg *.log *.out 
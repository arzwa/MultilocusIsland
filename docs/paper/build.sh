pandoc -o draft.pdf draft.md --bibliography /home/arthur_z/vimwiki/bib.bib -F pandoc-crossref -C -N --csl=genetics.csl --template template.latex
pandoc -o draft.tex draft.md --bibliography /home/arthur_z/vimwiki/bib.bib -F pandoc-crossref -C -N --csl=genetics.csl --template template.latex
#pandoc -o draft.pdf draft.md --bibliography /home/arthur_z/vimwiki/bib.bib -F pandoc-crossref -C -N --csl=plos-genetics.csl --template template.latex


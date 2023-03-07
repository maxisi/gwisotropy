import paths
from parse import parse
from glob import glob
import os
import numpy as np

skypath = str(paths.figures / "skymap_j_{}.pdf")
fname = os.path.basename(skypath)

events = []
for p in sorted(glob(skypath.format('*'))):
    events.append(parse(skypath, p)[0])

ncol = 5
nrow = 9

ntot = ncol*nrow
width = 1./ncol

npage = int(np.ceil(len(events) / ntot))

s = r"""\begin{subfigure}{%.2f\textwidth}
    \centering
    \includegraphics[width=\textwidth]{figures/%s}
    \caption{%s}
\end{subfigure}"""

header = r"""\begin{figure*}
\captionsetup[subfigure]{labelformat=empty}
"""

footer = r"""
\caption{Measurements of the total angular momentum direction, $\hat{J}$, for
the events in our set, in a Molweide projection of Earth-centric Celestial coordinates; darker color represents higher probability density for that direction in space%s.}
\label{fig:skymaps-%i}
\end{figure*}"""

outpath = str(paths.output / "skymaps.tex")

figures = []
for ipage in range(npage):
    figure = header
    for i, e in enumerate(events[ipage*ntot:(ipage+1)*ntot]):
        subfig = s % (width, fname.format(e), e.replace('_', '\\_'))
        if (i + 1) % ncol:
            subfig += "%"
        figure += subfig + '\n'
    figure = figure.strip('%\n')
    if ipage > 0:
        # add ContinuedFloat to make sure all figures get same Fig. number 
        figure = figure.replace(r"\begin{figure*}",
                                r"\begin{figure*}\ContinuedFloat")
        # add (continued) to caption
        note = " (cont.)"
    else:
        note = ""
    # add footer to figure
    figure += footer % (note, ipage+1)
    if npage > 1 and ipage < npage-1:
        figure += "%\n"
    figures.append(figure)
    
with open(outpath.format(ipage), 'w') as f:
    f.write(''.join(figures))

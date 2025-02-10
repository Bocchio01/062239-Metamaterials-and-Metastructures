function export_pdf_graphic(tile, path)

if (nargin == 1)
    path = '/untitled';
end

filename = ['.latex/report/img/MATLAB' path '.pdf'];
exportgraphics(tile, filename, 'ContentType', 'vector');

end


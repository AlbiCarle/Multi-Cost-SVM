function cool_image_to_pdf(inputnamefig, outputnamefig)

% Load the .fig file
figFileName = strcat(inputnamefig,'.fig');  
hFig = openfig(figFileName);

% Get the handle to the axis in the figure
hAxis = gca;

% Set the axis to have tight margins
set(hAxis, 'LooseInset', get(hAxis, 'TightInset')+40);

% Save the figure as a cropped PDF without any additional margins
outputFileName = strcat(outputnamefig,'.pdf');   
exportgraphics(hFig, outputFileName, 'ContentType', 'vector');

% Close the figure handle
close(hFig);

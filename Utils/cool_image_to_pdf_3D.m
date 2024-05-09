function cool_image_to_pdf_3D(inputnamefig, outputnamefig)

% Load the .fig file
figFileName = strcat(inputnamefig, '.fig');
hFig = openfig(figFileName);

% Get the handle to the axis in the figure
hAxis = gca;

% Set the axis to have tight margins
set(hAxis, 'LooseInset', get(hAxis, 'TightInset')+0.01);

% Define the output image file name (e.g., PNG)
outputImageFileName = strcat(outputnamefig, '.png');

% Export the figure as a PNG image
print(hFig, outputImageFileName, '-dpng', '-r300'); % Adjust resolution as needed

% Close the figure handle
close(hFig);
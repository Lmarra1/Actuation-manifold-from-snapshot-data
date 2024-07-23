function ReadH5(filename)

    % Retrieve information about the HDF5 file
    info = h5info(filename);

    % Loop through each dataset in the file and read the data
    for i = 1:length(info.Datasets)
        datasetName = info.Datasets(i).Name;
        variableName = matlab.lang.makeValidName(datasetName);  % Generate a valid MATLAB variable name
        
        % Read the data from the dataset
        data = h5read(filename, ['/' datasetName]);
        
        % Dynamically create a variable in the base workspace
        assignin('base', variableName, data);
    end
end



% info = h5info(filename);
% 
% % Initialize a structure to store the data
% data = struct();
% 
% % Loop through each dataset in the file and read the data
% for i = 1:length(info.Datasets)
%     datasetName = info.Datasets(i).Name;
%     variableName = genvarname(datasetName);  % Generate a valid MATLAB variable name
%     
%     % Read the data and assign it to a field with the generated name
%     data.(variableName) = h5read(filename, ['/' datasetName]);
% end

% end
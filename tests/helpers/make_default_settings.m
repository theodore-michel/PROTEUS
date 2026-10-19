function S = make_default_settings(PATHS, geometryFolder)
%MAKE_DEFAULT_SETTINGS Headless equivalent of the MainGUI start-up.
%
% S = MAKE_DEFAULT_SETTINGS(PATHS) returns a struct with the fields
% Microbubble, Medium, Acquisition, Transducer, Transmit,
% SimulationParameters and Geometry, initialised exactly as MainGUI does
% in startupFcn followed by updateParameters.
%
% S = MAKE_DEFAULT_SETTINGS(PATHS, geometryFolder) selects another flow
% geometry folder in PATHS.GeometriesPath (default: 'renal_tree').

if nargin < 2
    geometryFolder = 'renal_tree';
end

S.Microbubble          = reset_microbubble();
S.Medium               = reset_medium();
S.Acquisition          = reset_acquisition();
S.Transducer           = reset_transducer();
S.Transmit             = reset_transmit(S.Transducer);
S.SimulationParameters = reset_simulation_parameters();

Geometry.Folder = geometryFolder;
S.Geometry = reset_geometry(Geometry, PATHS);

S = update_derived_settings(S);
end

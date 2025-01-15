% clc
% clear variables
% close all

warning('off', 'MATLAB:print:ContentTypeImageSuggested');

%% Simulation name decode

modulation_model = [
    "ON-ON-ON"
    "OFF-OFF-OFF"
    "ON-ON-OFF"
    ]';

Z_su_model = [
    "Absent" 
    "R"
    "L"
    "C-"
    "RL"
    "RC-"
    "RLC-"
    ]';

E_su_model = [
    "Generic"
    "C- (Ideal)"
    "C- (Real)"
    ]';

simulation_type = [
    "TM"
    "PWEM"
    "FDTD"
    ]';


%% App launcher

app(modulation_model, Z_su_model, E_su_model, simulation_type);
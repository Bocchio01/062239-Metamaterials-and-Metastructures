% clc
% clear variables
% close all

warning('off', 'MATLAB:print:ContentTypeImageSuggested');

%% Simulation name decode

STC_modulation_label = [
    "ON-ON-ON"
    "OFF-OFF-OFF"
    "ON-ON-OFF"
    "Sinusoidal (continuos)"
    "Sinusoidal (discrete)"
    ]';

Z_model_label = [
    "Absent" 
    "R"
    "L"
    "C-"
    "RL"
    "RC-"
    "RLC-"
    ]';

simulation_label = [
    "TM"
    "PWEM"
    "FDTD"
    ]';


%% App launcher

app(STC_modulation_label, Z_model_label, simulation_label);
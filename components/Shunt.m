classdef Shunt
    % SHUNT Class representing a shunt impedance circuit.
    % This class models an electrical shunt circuit with resistive,
    % inductive, and capacitive elements, providing an impedance function
    % based on a given configuration label.

    properties
        R0  % High-value resistance (Ohms)
        R1  % Series resistance 1 (Ohms)
        R2  % Series resistance 2 (Ohms)
        C0  % Capacitance C0 (Farads)
        L   % Inductance (Henries)
        R   % Resistance (Ohms)
        C   % Additional capacitance (Farads)
        Z   % Function handle for impedance calculation
    end

    properties (SetAccess = immutable)
        label  % Label specifying the shunt configuration (e.g., 'RLC-', 'RC//L')
    end

    methods
        function obj = Shunt(label, varargin)
            % Constructor for the Shunt class.
            % Initializes the circuit parameters and computes impedance.
            %
            % Syntax:
            %   obj = Shunt(label, 'R0', value, 'R1', value, ...)
            %
            % Inputs:
            %   label - (string) Specifies the shunt configuration.
            %   Name-Value pairs:
            %       'R0' - Resistance R0 (default: 1000 kΩ)
            %       'R1' - Resistance R1 (default: 7.5 kΩ)
            %       'R2' - Resistance R2 (default: 13.7 kΩ)
            %       'C0' - Capacitance C0 (default: 4.4 nF)
            %       'L'  - Inductance (default: 15 mH)
            %       'R'  - Resistance R (default: 1000 Ω)
            %       'C'  - Additional capacitance (default: 5 nF)
            %
            % Outputs:
            %   obj - Instance of the Shunt class.

            obj.label = label; % Store the configuration label

            % Create an input parser for optional parameters
            parser = inputParser;
            addParameter(parser, 'R0', 1000e3);  % Default 1 MΩ
            addParameter(parser, 'R1', 7.5e3);   % Default 7.5 kΩ
            addParameter(parser, 'R2', 13.7e3);  % Default 13.7 kΩ
            addParameter(parser, 'C0', 4.4e-9);  % Default 4.4 nF
            addParameter(parser, 'L', 15e-3);    % Default 15 mH
            addParameter(parser, 'R', 1000);     % Default 1000 Ω
            addParameter(parser, 'C', 5e-9);     % Default 5 nF

            % Parse the inputs
            parse(parser, varargin{:});

            % Assign parsed values to object properties
            obj.R0 = parser.Results.R0;
            obj.R1 = parser.Results.R1;
            obj.R2 = parser.Results.R2;
            obj.C0 = parser.Results.C0;
            obj.L = parser.Results.L;
            obj.R = parser.Results.R;
            obj.C = parser.Results.C;

            % Compute the impedance function based on the given label
            obj.Z = obj.getZModel();
        end

        function Z = getZModel(obj)
            % GETZMODEL Returns the impedance function for the given shunt configuration.
            %
            % This function maps different shunt configurations (label) to their
            % corresponding impedance function handle.
            %
            % Outputs:
            %   Z - Function handle @(w) Z(w) representing impedance as a function of frequency.

            % Define a map linking circuit configurations to impedance functions
            Z_map = containers.Map( ...
                {'OFF', '+inf', 'C-', 'RLC-', 'RLC', 'RL//C', 'RC//L'}, ...
                { ...
                    @(w) 0, ...  % OFF: Zero impedance
                    @(w) +inf, ...  % Infinite impedance
                    @(w) - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w * obj.C0)^-1, ... % C-
                    @(w) obj.R + 1i*w * obj.L - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w * obj.C0)^-1, ... % RLC-
                    @(w) obj.R + 1i*w * obj.L + (1i*w * obj.C)^-1, ... % RLC
                    @(w) (obj.R + 1i*w * obj.L) / (1 + 1i*w * obj.R*obj.C - w^2 * obj.L*obj.C), ... % RL//C
                    @(w) (1i*w * obj.R*obj.L + obj.L / obj.C) / (1 + 1i*w * obj.R*obj.C - w^2 * obj.L*obj.C) ... % RC//L
                } ...
                );

            % Retrieve the impedance function for the given label
            Z = Z_map(obj.label);
        end
    end
end

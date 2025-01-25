classdef Shunt

    properties
        R0
        R1
        R2
        C0
        L
        R
        Z
    end

    properties (SetAccess = immutable)
        label
    end
    
    methods
        function obj = Shunt(label, varargin)

            obj.label = label;

            parser = inputParser;
            addParameter(parser, 'R0', 1000e3);
            addParameter(parser, 'R1', 7.5e3);
            addParameter(parser, 'R2', 13.7e3);
            addParameter(parser, 'C0', 4.4e-9);
            addParameter(parser, 'L', 15e-3);
            addParameter(parser, 'R', 1000);
            
            % Parse the inputs
            parse(parser, varargin{:});
            
            % Assign properties
            obj.R0 = parser.Results.R0;
            obj.R1 = parser.Results.R1;
            obj.R2 = parser.Results.R2;
            obj.C0 = parser.Results.C0;
            obj.L = parser.Results.L;
            obj.R = parser.Results.R;

            obj.Z = obj.getZModel();

        end
        
        function Z = getZModel(obj)
            % We need to complete with the paraller connection of the shunt
            % elements R L C-
            Z_map = containers.Map( ...
            {'Absent', 'R', 'L', 'C-', 'RL', 'RC-', 'RLC-'}, ...
            { ...
                @(w) +inf, ...
                @(w) obj.R, ...
                @(w) 1i*w*obj.L, ...
                @(w) - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w*obj.C0)^-1, ...
                @(w) obj.R + 1i * obj.L * w, ...
                @(w) obj.R - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w*obj.C0)^-1, ...
                @(w) obj.R + 1i * obj.L * w - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w*obj.C0)^-1 ...
            } ...
            );

            Z = Z_map(obj.label);

        end
    end
end


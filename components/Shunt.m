classdef Shunt

    properties
        R0
        R1
        R2
        C0
        L
        R
        C
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
            addParameter(parser, 'C', 5e-9);
            
            % Parse the inputs
            parse(parser, varargin{:});
            
            % Assign properties
            obj.R0 = parser.Results.R0;
            obj.R1 = parser.Results.R1;
            obj.R2 = parser.Results.R2;
            obj.C0 = parser.Results.C0;
            obj.L = parser.Results.L;
            obj.R = parser.Results.R;
            obj.C = parser.Results.C;

            obj.Z = obj.getZModel();

        end
        
        function Z = getZModel(obj)
            % We need to complete with the parallel connection of the shunt
            % elements R L C-
            Z_map = containers.Map( ...
            {'OFF', '+inf', 'C-', 'RLC-', 'RLC', 'RL//C', 'RC//L'}, ...
            { ...
                @(w) 0, ...
                @(w) +inf, ...
                @(w) - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w * obj.C0)^-1, ... % C-
                @(w) obj.R + 1i*w * obj.L - obj.R1 / obj.R2 * (1/obj.R0 + 1i*w * obj.C0)^-1, ... % RLC-
                @(w) obj.R + 1i*w * obj.L + (1i*w * obj.C)^-1 ... % RLC
                @(w) (obj.R + 1i*w * obj.L) / (1 + 1i*w * obj.R*obj.C - w^2 * obj.L*obj.C) ...
                @(w) (1i*w * obj.R*obj.L + obj.L / obj.C) / (1 + 1i*w * obj.R*obj.C - w^2 * obj.L*obj.C) ...
            } ...
            );

            Z = Z_map(obj.label);

        end
    end
end


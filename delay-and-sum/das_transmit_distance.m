function d_in = das_transmit_distance(X, Z, transducer, transmit)
%DAS_TRANSMIT_DISTANCE Transmit path length of the delay-and-sum model.
%
% d_in = DAS_TRANSMIT_DISTANCE(X, Z, transducer, transmit) returns, for
% each point (X, Z) (lateral, axial) [m], the distance travelled by the
% transmitted wave from the moment the first element fires (t = 0, the
% smallest transmit delay) until the wave reaches the point. d_in has the
% size of X.
%
% transducer: struct with fields
% - x_el, z_el:  lateral and axial element positions [m]
% - center:      axial coordinate of the centre of curvature (-Inf for a
%                linear array) [m]
% transmit: struct with fields
% - focus:  focal distance [m]: > 0 focused, < 0 diverging (virtual
%           source behind the array), Inf: plane wave for a linear array,
%           diverging wave from the centre of curvature for a curved array
% - steer:  steering angle [rad], positive towards positive lateral
%           coordinates
%
% The model matches the PROTEUS transmit delays (compute_delays): for each
% element, d_in at its own position equals c times its delay. As in
% compute_delays, t = 0 is the smallest delay over all elements, whether
% they fire (apodization) or not.

x_el = transducer.x_el(:)';
z_el = transducer.z_el(:)';

F     = transmit.focus;
steer = transmit.steer;

if ~isfinite(F) && isfinite(transducer.center)
    % Curved array: diverging wave from the centre of curvature
    F = transducer.center;
end

if isfinite(F)
    % Focal point or virtual source:
    Fx = F*sin(steer);
    Fz = F*cos(steer);
    D_el = sqrt((x_el-Fx).^2 + (z_el-Fz).^2);
    D    = sqrt((X-Fx).^2 + (Z-Fz).^2);

    if F < 0
        % Diverging wave, t = 0 at the element closest to the virtual
        % source:
        d_in = D - min(D_el);
    else
        % Converging wave, t = 0 at the element farthest from the
        % focus. Points before the focus (s = 1) are reached earlier than
        % the focus, points after it (s = -1) later:
        s = sign((Fx-X)*sin(steer) + (Fz-Z)*cos(steer));
        d_in = max(D_el) - s.*D;
    end
else
    % Plane wave, t = 0 at the first element along the propagation
    % direction:
    d_in = X*sin(steer) + Z*cos(steer) - ...
        min(x_el*sin(steer) + z_el*cos(steer));
end

end

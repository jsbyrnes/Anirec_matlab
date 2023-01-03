function [ P_comp, SV_comp, SH_comp ] = anirec( phase, dt, slow, baz, model )
%ANIREC Runs the anirec program, and returns a normalized three component green's
%function in ray coordinates. Length of returned trace is fixed to 100 seconds. Phase is case
%insensitive. dt in seconds, slow in sec/km, baz in degree clockwise from
%north. Both velocities in km/s, z is depth to base of layer. Last layers
%is infinite half space. 

    phase = upper(phase);

    if ~strcmpi(phase, 'P') && ~strcmpi(phase, 'SV')
        
        error('Phase should be P or SV');
        
    end
    
    model.vp = model.vp*1000;
    model.vs = model.vs*1000;
    model.z = model.z*1000;
    model.rho = model.rho*1000;
    
    cc = 1/slow;
    
    if strcmpi(phase, 'P')        

        [z_comp, r_comp, SH_comp] = anirec_gateway(model.theta, model.phi, model.z, model.vp, ...
            model.A, model.B, model.vs, model.C, model.rho, length(model.vp) - 1, cc, baz, dt, 1);
        
        [P_comp, SV_comp] = ZR2PSV(z_comp, r_comp, slow, model.vp(1)/1000);
        
    elseif strcmpi(phase, 'SV')
        
        [z_comp, r_comp, SH_comp] = anirec_gateway(model.theta, model.phi, model.z, model.vp, ...
            model.A, model.B, model.vs, model.C, model.rho, length(model.vp) - 1, cc, baz, dt, 2);
        
        SV_comp = -1*r_comp;
        P_comp = -1*z_comp;
        
        %[SV_comp, P_comp] = ZR2PSV(r_comp, z_comp, slow, model.vs(1)/1000);
        
    end
    
end


function [ source_comp_after conv_comp_after ] = ZR2PSV( source_comp_before, conv_comp_before, rp, surfv)
%ZR2PSV Rotate the Z and R traces into ray coordinates
%   surfv should be the vp(1) for a P and vs(1) for an S
%   Make sure that rp is plane wave(s/km). Source component is the 
%   channel that the incoming wave in on, conv component is the channe;
%   that the converted phases are on. Also renormalizes after rotation
%   to incoming component. Assumes all later phases are less than incoming
%   pulse(probably true for synthetic data, want to be true for real data)
%
%   Joseph Byrnes
%   Oct 15th 2012, jbyrnes@uoregon.edu
%

    angle = asin(surfv*rp);

    conv_comp_after = cos(angle)*conv_comp_before - sin(angle)*source_comp_before;
    source_comp_after = sin(angle)*conv_comp_before + cos(angle)*source_comp_before;
   
    amp = max(source_comp_after);
    
    source_comp_after = source_comp_after/amp;
    conv_comp_after = conv_comp_after/amp;
     
end
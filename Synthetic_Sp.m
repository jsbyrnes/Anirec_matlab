%This script make a synthetic receiver function filtered at a given band.
%That's it.

clear, close all

addpath('./toolbox/');
addpath('./deconvolution_code');

filter_ends = [ 0.01 0.25];
source      = [ 0.25 0.2 ];

figure(1)

dz   = 1;
vpvs = 1.76;

fid = fopen('./Data7621_bw1_05025.40.0.-103.0.csv');
fgetl(fid);fgetl(fid);fgetl(fid);
bids = fgetl(fid);
fclose(fid);

bids = strsplit(bids, ',');

indM = str2num(bids{1}) + 2;
indN = str2num(bids{2}) + 2;

m = readmatrix('./Data7621_bw1_05025.40.0.-103.0.csv', 'NumHeaderLines', 4);

m(:, 2) = cumsum(m(:, 2));

Ndepth = 0.5*(m(indN, 2) + m(indN-1, 2));

subplot(1,2,1)
plot(m(:, 1), m(:, 2), 'k', 'LineWidth', 2)
xlabel('Vs, km/s'), ylabel('Depth, km')
set(gca, 'YDir', 'reverse')
ylim([20 120])
grid on

%%%%%%%%%%%
%make a new model in the format of the model structure
depth = dz:dz:200;

model.z    = depth;
model.vp   = vpvs*interp1(m(:, 2), m(:, 1), depth);
model.vpvs = vpvs*ones(size(model.vp));

%now fill out with zeros, not doing anisotropy here, then density
model.theta = zeros(size(model.vp));
model.phi   = zeros(size(model.vp));
model.A     = zeros(size(model.vp));
model.B     = zeros(size(model.vp));
model.C     = zeros(size(model.vp));
model.rho   = nafedrake_rho(model.vp); %g/cm^3

%now set other parameters needed to make a receiver function
baz    = 0;
shift  = 40;
phase  = 'SV';
slow   = .05;
dt     = .05;

[P_comp, SV_comp, SH_comp] = anirec(phase, dt, slow, baz, model, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make a source that looks like a teleseismic S wave in frequency
Nfft                   = length(P_comp);         % number  of points in fft = n of samples
dF                     = 1/(dt*Nfft);       % frequency interval
Nyq                    = (1/(2*dt));        % Nyquist frequency
freqVals               = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);
freqVals               = abs(freqVals);

Faux = source(2)/sqrt(log(2));%number is half width in Hz
A    = exp(-((freqVals-source(1))/Faux).^2);
A    = A'; %for dimensional consistency with inTr.data

%Apply to the traces and revert to time domain
P_comp  = real(ifft(fft(P_comp).*A));
SV_comp = real(ifft(fft(SV_comp).*A));
SH_comp = real(ifft(fft(SH_comp).*A));

[trace_phase, ~] = multitaper2rf_3component(P_comp, SV_comp, SH_comp, dt, shift, 2000, 2.5, 3, phase, filter_ends, [ 1 length(SV_comp)], []);
t = ( 0:dt: (length(trace_phase) - 1)*dt) - shift;

%migrate to depth by finding the S-p time through each layer
zt = zeros(size(depth));
for k = 1:length(zt)

    %s - p time for actual receiver function
    zt(k)    = dz*(sqrt(1/(model.vp(k)/vpvs)^2 - slow^2) - sqrt(1/(model.vp(k))^2 - slow^2));
    %vertical s time for comparing to Emily's data

    if k < length(zt)

        stime(k) = dz/(0.5*model.vp(k) + 0.5*model.vp(k+1));

    else

        stime(k) = 0;%not correct but whatever

    end

end
zt    = cumsum(zt);
stime = cumsum(stime)*vpvs;

trace_phase = interp1(t, trace_phase, zt);
%now simulate the pick
%[pks,locs,widths,prom] = findpeaks(-1*trace_phase,depth, 'WidthReference', 'halfheight');
[pks,locs,widths,prom] = findpeaks(-1*trace_phase,depth);

[~, ind] = min( abs(locs - Ndepth));
%[~, ind] = max(pks);

%cull to the peak, then cull to halfwidth

w = trace_phase+pks(ind);

pkinds = (depth < (locs(ind) + widths(ind)/2)) & (depth > (locs(ind) - widths(ind)/2)) & (w < prom(ind)/2);

depth_pick = sum(w(pkinds).*depth(pkinds))/sum(w(pkinds));

predicted_time = interp1(depth, stime, depth_pick);

subplot(1,2,2)
hold on
plot(trace_phase, depth, 'Linewidth', 2)
plot([-1 1], [depth_pick depth_pick], 'k--');
xlabel('RF amplitude')
set(gca, 'YDir', 'reverse')
xlim( [(min(trace_phase) - 0.01) (max(trace_phase) + 0.01)]);
grid on
ylim([20 120])

%and now put it on the velocity model
subplot(1,2,1)
hold on
plot([-10 10], [depth_pick depth_pick], 'k--');
xlim( [(min(m(:, 1)) - 0.1) (max(m(:, 1)) + 0.1)]);

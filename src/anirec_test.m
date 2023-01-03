clear, close all
model.theta = [ 90 0 0];%tilt of the fast azis, zero is verticle
model.phi = [ 0 0 0 ];%azimuth of the fast direction
model.z = [30 90 300];
model.vp = [8 8 8];
model.A = [ 0.1 0 0 ];
model.B = [0.0 0 0 ];
model.vs = model.vp/1.76;
model.C = [ 0.1 0 0 ];
model.rho = .33 + .77*model.vp;

phase = 'SV';

baz = 45;

slow = .04;

dt = .05;

tic

[Pcomp, SVcomp, SHcomp] = anirec_shtest(phase, dt, slow, baz, model);

toc

%Pcomp = [ Pcomp;zeros(4500-2001,1)];
%SVcomp = [ SVcomp;zeros(4500-2001,1)];
%SHcomp = [ SHcomp;zeros(4500-2001,1)];

% Pcomp = circshift(Pcomp, 2300);
% SVcomp = circshift(SVcomp, 2300);
% SHcomp = circshift(SHcomp, 2300);

% figure, plot(real(Pcomp)), figure, plot(real(SVcomp), 'k'), figure, plot(real(SHcomp), 'b')
% 
% figure, hold on
% plot(real(Pcomp), 'k');
% plot(real(SVcomp), 'b');
% plot(imag(hilbert(SVcomp)), 'b--');

figure(1)
plot(Pcomp)
title('P component');

figure(2)
plot(SVcomp);
title('SV component');

figure(3)
plot(SHcomp)
title('SH component');
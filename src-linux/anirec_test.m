clear, close all
model.theta = [90 0 0];
model.phi = [0 0 0 ];
model.z = [25 90 300];
model.vp = [8 8 8];
model.A = [ 0 0 0 ];
model.B = [0 0 0 ];
model.vs = model.vp/1.76;
model.C = [ 0 0 0 ];
model.rho = .33 + .77*model.vp;

phase = 'P';

baz = 45;

slow = .2;

dt = .05;

tic

[Pcomp, SVcomp, SHcomp] = anirec(phase, dt, slow, baz, model);

toc

%Pcomp = [ Pcomp;zeros(4500-2001,1)];
%SVcomp = [ SVcomp;zeros(4500-2001,1)];
%SHcomp = [ SHcomp;zeros(4500-2001,1)];

% Pcomp = circshift(Pcomp, 2300);
% SVcomp = circshift(SVcomp, 2300);
% SHcomp = circshift(SHcomp, 2300);

figure, plot(real(Pcomp)), figure, plot(real(SVcomp), 'k'), figure, plot(real(SHcomp), 'b')

figure, hold on
plot(real(Pcomp), 'k');
plot(real(SVcomp), 'b');
plot(imag(hilbert(SVcomp)), 'b--');
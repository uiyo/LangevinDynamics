clear all
close all

% Driver script
% author: Yuguang Yang yyang60@jhu.edu

% first setting simulation parameter
SetSimPara;
rng(1);
% run simulation
x = Simulator(simPara);
time = [1:size(x,2)]*simPara.dt*simPara.saveInterval;

% plot the trajectory
figure(1)
plot(x(1,:)/simPara.radius,x(2,:)/simPara.radius)
xlabel('x/a')
ylabel('y/a')

figure(2)
plot(time,x(3,:)/simPara.radius)
xlabel('t,s')
ylabel('z/a')
figure(3)
plot(time,x(1,:)/simPara.radius)
xlabel('t,s')
ylabel('x/a')
figure(4)
plot(time,x(2,:)/simPara.radius)
xlabel('t,s')
ylabel('y/a')
figure(5)
hist(x(3,:)/simPara.radius,100)
xlabel('z/a')
ylabel('count')
[n c] = hist(x(3,:)/simPara.radius,100);
figure(6)
u = -log(n/sum(n));
u = u - min(u);
plot(c,u);
ylim([0 7])
xlabel('z/a')
ylabel('u/kT')

mvflag = 1;

if mvflag == 1
    len = length(x);
    skip = 2;
    for i=1:skip:len/5
        figure(10)
        plot(x(1,i)/simPara.radius,x(2,i)/simPara.radius,'linestyle','none','marker','o','markersize',8);
         xlim([-0.2 0.2])
         ylim([-1 1])
         saveas(gcf,['movie3/frame' num2str(i,'%03d') '.jpg'])
        
    end
end

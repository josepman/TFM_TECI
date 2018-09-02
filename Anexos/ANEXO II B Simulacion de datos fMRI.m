
sP = simtb_create_sP('experiment_params_block_mod');

simtb_main(sP);

openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_model.fig')

datos = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_001_SIM.mat')

openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_001.fig')

openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_SM_translate_x.fig')
openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_SM_translate_y.fig')
openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_SM_theta.fig')
openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_SM_spread.fig')

f = openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/block_D_pSC.fig', figsize=(10,10))

A = imread('/Users/hose/Desktop/TFM_TECI/simulated_data/aod_motion005.jpg');
image(A)

simtb_showTC(sP);

off = zeros(1,30);
on = ones(1,20);
task = [off,on,off,on,off,on,off,on,off,on];
%csvwrite('task_design.csv', task)

plot(task, 'k--' ), hold on, 
plot(TC_block(:,16), 'g', 'linewidth',2), hold on, 
plot(TC_block(:,22), 'r', 'linewidth',2), hold on, 
plot(TC_block(:,23), 'm', 'linewidth',1), hold on, 
plot(TC_block(:,6),'b','linewidth',2)
axis([0 250 -1 1.3])
legend('Task', 'Motor', 'Auditivo', 'Hipocampo', 'DMN')
xlabel('Tiempo (en TRs)')
ylabel('Amplitud')
title('Tarea por bloques y respuesta')

openfig('/Users/hose/Desktop/TFM_TECI/simulated_data/task.fig')

s = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_001_SIM.mat');
tc1 = s.TC;
sm1 = s.SM;
cmtc1 = s.cmTC;
cmsm1 = s.cmSM;
s = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_002_SIM.mat');
tc2 = s.TC;
sm2 = s.SM;
cmtc2 = s.cmTC;
cmsm2 = s.cmSM;
s = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_003_SIM.mat');
tc3 = s.TC;
sm3 = s.SM;
cmtc3 = s.cmTC;
cmsm3 = s.cmSM;
s = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_004_SIM.mat');
tc4 = s.TC;
sm4 = s.SM;
cmtc4 = s.cmTC;
cmsm4 = s.cmSM;
s = load('/Users/hose/Desktop/TFM_TECI/simulated_data/block_subject_005_SIM.mat');
tc5 = s.TC;
sm5 = s.SM;
cmtc5 = s.cmTC;
cmsm5 = s.cmSM;

tc_avg = (tc1+tc2+tc3+tc4+tc5)/5;
sm_avg = (sm1+sm2+sm3+sm4+sm5)/5;
cmtc_avg = (cmtc1+cmtc2+cmtc3+cmtc4+cmtc5)/5;
cmsm_avg = (cmsm1+cmsm2+cmsm3+cmsm4+cmsm5)/5;

pcolor(cmtc_avg)

t = [1:200];
figure1=figure('Position', [1000, 1000, 1024, 1200]);
subplot(331), plot(t,tc_avg(t,3)), title('Frontal L(IC=4)') 
subplot(332), plot(t,tc_avg(t,4)), title('Frontal R(IC=4)')  
subplot(333), plot(t,tc_avg(t,5)), title('Precuneus')
subplot(334), plot(t,tc_avg(t,6)), title('DMN (IC=8)')    % Default Mode Network
subplot(335), plot(t,tc_avg(t,16)), title('Motor (IC=22)')
subplot(336), plot(t,tc_avg(t,17)), title('Motor (IC=23)')
subplot(337), plot(t,tc_avg(t,21)), title('Auditivo (IC=27)')
subplot(338), plot(t,tc_avg(t,22)), title('Auditivo (IC=28)')
subplot(339), plot(t,tc_avg(t,23)), title('Hipocampo (IC=29)')

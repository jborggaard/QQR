  setKroneckerSumPath

%  Set up test examples, problem dimensions (order), and degree of feedback

  addpath('./examples')
  addpath('./examples/ChafeeInfante1DControl/')
  addpath('./examples/KuramotoSivashinsky1DControl/')
  
%  Burgers equation


% %  Chafee-Infante equation
%  n = 24;
%  m = 10;
%  degree = 5;
%  a = 1e-3;    % scaling parameter on initial conditions
%  controlType = 'source';
%  r0=100.0;
%   
%  exampleCI
%  save('CIm10d5n41R100src_01.mat','k','v')
%  save('CIm10d5n41R100src_02.mat','t1','z1')
%  save('CIm10d5n41R100src_03.mat','t3','z3','a','controlType','r0')
%  save('CIm10d5n41R100src_04.mat','t5','z5')
 
%  Kuramoto-Sivashinsky equation
  n = 41;
  m = 5;
  degree = 4;

  L=13.5;
  tInfinity = 120;
  r0 = 1.0;
  exampleKS
  
  save('KSm5d4n41R001.mat')

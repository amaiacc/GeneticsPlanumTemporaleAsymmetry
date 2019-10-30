close all
clear all

load hospital;
hospital.SexID = grp2idx(hospital.Sex);
x = [hospital.SexID hospital.Age hospital.Smoker hospital.Weight];
rho = partialcorr(x);
rho = array2table(rho, ...
    'VariableNames',{'SexID','Age','Smoker','Weight'},...
    'RowNames',{'SexID','Age','Smoker','Weight'});

disp('Partial Correlation Coefficients')
disp(rho)


rho_test(1,1) = 1;
rho_test(2,2) = 1;
rho_test(3,3) = 1;
rho_test(4,4) = 1;
rho_test(1,2) = partialcorr(x(:,1),x(:,2),  [x(:,3), x(:,4)]);
rho_test(2,1) = rho_test(1,2);
rho_test(1,3) = partialcorr(x(:,1),x(:,3),  [x(:,2), x(:,4)]);
rho_test(3,1) = rho_test(1,3);
rho_test(1,4) = partialcorr(x(:,1),x(:,4),  [x(:,2), x(:,3)]);
rho_test(4,1) = rho_test(1,4);
rho_test(2,3) = partialcorr(x(:,2),x(:,3),  [x(:,1), x(:,4)]);
rho_test(3,2) = rho_test(2,3);
rho_test(2,4) = partialcorr(x(:,2),x(:,4),  [x(:,1), x(:,3)]);
rho_test(4,2) = rho_test(2,4);
rho_test(3,4) = partialcorr(x(:,3),x(:,4),  [x(:,1), x(:,2)]);
rho_test(4,3) = rho_test(3,4);